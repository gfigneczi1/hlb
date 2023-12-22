import numpy as np
import scipy.io
from scipy import interpolate
from _converter.main import Application
import os
from os import path
import logging
import read_rosbag
from glob import glob
from operator import itemgetter
import json
from math import ceil


class Convert:

    def __init__(self, extension, measurement_files, root, config_signals, vehicle_profile):
        self.extension = extension
        self.measurement_files = measurement_files
        self.root = root
        self.config_signals = config_signals
        self.vehicle_profile = vehicle_profile

    def convert_handler(self):
        if self.extension != "bag" and self.extension != "mat":
            self.convert_mdf()
        elif self.extension == "bag":
            self.convert_bag()
        elif self.extension == "mat":
            self.convert_mat()

    def convert_mdf(self):
        Application(self.vehicle_profile)

    def convert_bag(self):
        logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
        temp_folder_path = path.join(self.root, '_temp')
        extracted_files = []

        rosbag_reader = read_rosbag.initialize()
        file_ID = 0
        for meas_file in self.measurement_files:
            try:
                if file_ID > 0:
                    with open(path.join(self.root, '_temp', 'config.json')) as f:
                        config_json = json.load(f)
                    self.config_signals = config_json["signal_list"]
                new_signal_names = list(map(itemgetter(0), self.config_signals))
                header_dict_keys = tuple(map(lambda x: x + "_header", [name[0] for name in new_signal_names]))
                signal_headers = tuple(map(lambda x: x[0], new_signal_names))
                data_field = list(map(itemgetter(1), self.config_signals))
                filename = meas_file.split(sep="\\")[-1]
                logging.info(f"Started extracting data from {filename}")
                ros_topics = rosbag_reader.read_rosbag(meas_file, temp_folder_path, nargout=1)
                ros_topics_msgs = dict.fromkeys(tuple(ros_topics))
                msgs_path = sorted(glob(path.join(temp_folder_path, '*.mat')), key=len)

                # loading saved topics to a dict
                index = 0
                for topic_file in msgs_path:
                    if topic_file.split(sep="\\")[-1] not in extracted_files:
                        topic_data = scipy.io.loadmat(topic_file)
                        ros_topics_msgs[ros_topics[index]] = topic_data['msgCellArray']
                        os.remove(topic_file)
                        index += 1
                logging.info("ROS topics are loaded, starting signal renaming with data extraction.")

                # creating a list of indexing depth of data
                # and creating keys for signal dict
                index_list = []
                indexes_to_remove = []
                new_signal_names_dimension = []
                k = 0
                for signal in data_field:
                    invalid_field_path = False
                    new_signal_names_dimension.append(1)
                    field_list = []
                    i = 0
                    while True:
                        data_for_dtype = ros_topics_msgs
                        if len(field_list) != 0:
                            field_copy = field_list
                            dim_list = []
                            for m in range(ros_topics_msgs[signal[0]].shape[0]):
                                data_for_dtype = ros_topics_msgs
                                for j in range(len(field_list)):
                                    if isinstance(field_list[j], str) and field_list[j][0] != "/":
                                        dim_list.append(data_for_dtype[field_list[j]].shape[1])
                                    if j == 1:
                                        data_for_dtype = data_for_dtype[m]
                                    else:
                                        data_for_dtype = data_for_dtype[field_list[j]]
                                data_dtype = data_for_dtype.dtype
                            if len(dim_list) > 0:
                                new_signal_names_dimension[k] = max(dim_list)
                        if len(field_list) == 0:
                            field_list.append(signal[i])
                            i += 1
                        elif len(data_dtype) != 0 and signal[i] in data_dtype.names:
                            field_list.append(signal[i])
                            i += 1
                        elif data_dtype.name == 'object' and i != len(signal):
                            field_list.append(0)
                        elif 'data_dtype' in locals() and len(data_dtype) == 0 and i == len(signal):
                            break
                        else:
                            logging.error("The given signal profile: " + str(signal) + " is invalid!")
                            invalid_field_path = True
                            indexes_to_remove.append(data_field.index(signal))
                            break
                    if not invalid_field_path:
                        index_list.append(field_list)
                    k += 1
                for index_to_remove in indexes_to_remove:
                    new_signal_names.pop(index_to_remove)
                    signal = data_field.pop(index_to_remove)[0]
                    ros_topics_msgs.pop(signal)
                if len([x for x in new_signal_names_dimension if x > 1]) > 0:
                    for i in range(len(new_signal_names)):
                        default_signal_name = new_signal_names[i][0]
                        if new_signal_names_dimension[i] > 1:
                            for j in range(new_signal_names_dimension[i]):
                                if j == 0:
                                    new_signal_names[i][j] = default_signal_name + "_" + str(j + 1)
                                else:
                                    new_name = default_signal_name + "_" + str(j + 1)
                                    new_signal_names[i].insert(j, new_name)
                            new_signal_names[i].insert(j + 1, (default_signal_name + "_" + "size"))
                new_signal_names = [item for sublist in new_signal_names for item in sublist]

                # extracting data to the new signal dict
                signal_dict = dict.fromkeys(new_signal_names)
                for i in range(len(data_field)):
                    size_name = list(map(itemgetter(0), self.config_signals))[i][-1]
                    if len(list(map(itemgetter(0), self.config_signals))[i]) == 1:
                        signal_array_size = len(list(map(itemgetter(0), self.config_signals))[i])
                    else:
                        signal_array_size = len(list(map(itemgetter(0), self.config_signals))[i]) - 1
                    if i == 0:
                        shift = 0
                    else:
                        shift += len(list(map(itemgetter(0), self.config_signals))[previous_i])
                    for j in range(len(ros_topics_msgs[data_field[i][0]])):
                        data_extracted = []
                        index_list[i][1] = j
                        data = ros_topics_msgs
                        for k in index_list[i]:
                            if isinstance(k, str) and k[0] != "/" and len(data.shape) > 0 and data.shape[1] > 1:
                                data_tree = data[k]
                                for m in range(data_tree.shape[1]):
                                    data = data_tree.item(m)
                                    if index_list[i].index(k) != len(index_list[i]) - 1:
                                        for n in index_list[i][index_list[i].index(k) + 1:]:
                                            data = data[n]
                                    data_extracted.append(data)
                            elif np.isscalar(data):
                                break
                            else:
                                data = data[k]
                        if len(data_extracted) == 0 and data.size > 0:
                            data_extracted.append(data)
                        if len(data_extracted) == 1:
                            actual_data_size = data_extracted[0].size
                        else:
                            actual_data_size = len(data_extracted)
                        if len(data_extracted) < signal_array_size:
                            for diff in range(signal_array_size - len(data_extracted)):
                                data_extracted.append(np.zeros((1, 1)))
                        for data_index in range(len(data_extracted)):
                            if data_extracted[data_index].ndim == 1:
                                data_extracted[data_index] = np.expand_dims(data_extracted[data_index], axis=1)
                            if signal_dict[new_signal_names[data_index + shift]] is None:
                                signal_dict[new_signal_names[data_index + shift]] = data_extracted[data_index]
                            else:
                                signal_dict[new_signal_names[data_index + shift]] = np.concatenate(
                                    (signal_dict[new_signal_names[data_index + shift]], data_extracted[data_index]), axis=1)
                            if data_index == len(data_extracted) - 1 and "size" in size_name:
                                if signal_dict[size_name] is None:
                                    signal_dict[size_name] = list()
                                signal_dict[size_name].append(actual_data_size)
                    previous_i = i
                for key in signal_dict.keys():
                    if not isinstance(signal_dict[key], list) and signal_dict[key].dtype.name == 'object':
                        signal_dict[key] = signal_dict[key].astype(np.float32)

                # extracting timestamps for signals
                logging.info("Interpolation started.")
                header_dict = dict.fromkeys(header_dict_keys)
                stamp_list = []
                for i in range(len(header_dict.keys())):
                    data_for_header = ros_topics_msgs
                    index_list[i][1] = 0
                    data_stamp = []
                    keys_to_delete = []
                    for j in index_list[i]:
                        data_for_header = data_for_header[j]
                        if hasattr(data_for_header,'dtype'):
                            if data_for_header.dtype.name is not 'object' and 'Header' in data_for_header.dtype.names:
                                data_stamp.extend([data_for_header['Header'][0][0]['Stamp'][0][0]['Sec'][0][0].item(),
                                            data_for_header['Header'][0][0]['Stamp'][0][0]['Nsec'][0][0].item()])
                                break
                    data_for_header = ros_topics_msgs
                    index_list[i][1] = -1
                    for j in index_list[i]:
                        data_for_header = data_for_header[j]
                        if hasattr(data_for_header,'dtype'):
                            if data_for_header.dtype.name is not 'object' and 'Header' in data_for_header.dtype.names:
                                data_stamp.extend([data_for_header['Header'][0][0]['Stamp'][0][0]['Sec'][0][0].item(),
                                                data_for_header['Header'][0][0]['Stamp'][0][0]['Nsec'][0][0].item()])
                                break
                    if len(data_stamp) == 4 and 0 not in data_stamp:
                        stamp_list.append(data_stamp)
                    else:
                        keys_to_delete = list(filter(lambda x: signal_headers[i] in x, signal_dict.keys()))
                        logging.warning("No valid timestamp found for data: {}, so it cannot be converted".format(keys_to_delete[:-1]))
                        for key in keys_to_delete:
                            signal_dict.pop(key)
                        header_dict.pop(header_dict_keys[i])
                if len(signal_dict) == 0:
                    logging.error("No valid data for interpolation, conversion aborted.")
                    break

                # interpolating the signals
                start_stamps = []
                end_stamps = []
                new_keys = list(header_dict.keys())
                for i in range(len(stamp_list)):
                    header_dict[new_keys[i]] = stamp_list[i]
                    start_stamps.append(float(".".join([str(stamp_list[i][0]), str(stamp_list[i][1])])))
                    end_stamps.append(float(".".join([str(stamp_list[i][2]), str(stamp_list[i][3])])))
                min_index, max_index = start_stamps.index(min(start_stamps)), end_stamps.index(max(end_stamps))
                q_T0_unix = [stamp_list[min_index][0], stamp_list[min_index][1],
                            stamp_list[max_index][2], stamp_list[max_index][3]]
                dt = 0.01
                q_T0_duration = float(end_stamps[max_index] - start_stamps[min_index])#float(".".join([str(q_T0_unix[2] - q_T0_unix[0]), str(q_T0_unix[3] - q_T0_unix[1])]))
                q_T0 = np.arange(0, q_T0_duration + dt, dt)
                for i in range(len(signal_dict)):
                    signal_dict_key = list(signal_dict.keys())[i]
                    if len(signal_dict_key) > 4 and "_size" in signal_dict_key:
                        pass
                    else:
                        old_data = signal_dict[signal_dict_key].flatten()
                        header_key = [key for key in list(header_dict.keys()) if key[:-7] in signal_dict_key][0]
                        data_duration = float(".".join([str(header_dict[header_key][2] - header_dict[header_key][0]),
                                                str(abs(header_dict[header_key][3] - header_dict[header_key][1]))]))
                        data_hz = round(data_duration / old_data.size, 3)
                        q_T0_unix_endtime = float(str(q_T0_unix[2])) + float(
                            str(q_T0_unix[3])) / 1000000000
                        q_T0_unix_starttime = float(str(q_T0_unix[0])) + float(
                            str(q_T0_unix[1])) / 1000000000
                        start_time = float(str(header_dict[header_key][0])) + float(str(header_dict[header_key][1])) / 1000000000 - q_T0_unix_starttime
                        end_time = q_T0_unix_endtime - (float(str(header_dict[header_key][2])) + float(
                            str(header_dict[header_key][3])) / 1000000000)

                        # start_sec = header_dict[header_key][0] - q_T0_unix[0]
                        # start_nsec = header_dict[header_key][1] - q_T0_unix[1]
                        # if start_nsec < 0:
                        #     start_sec -= 1
                        #     start_nsec = 1000000000 + start_nsec
                        # end_sec = q_T0_unix[2] - header_dict[header_key][2]
                        # end_nsec = q_T0_unix[3] - header_dict[header_key][3]
                        # if end_nsec < 0:
                        #     end_sec -= 1
                        #     end_nsec = 1000000000 + end_nsec

                        missing_start = ceil(start_time / data_hz)
                        missing_end = ceil(end_time / data_hz)
                        if missing_start < 0:
                            missing_start = 0
                        filled_data = np.insert(old_data, 0, np.zeros(missing_start))
                        filled_data = np.append(filled_data, np.zeros(missing_end), axis=None)

                        x = np.arange(0, q_T0_duration + data_hz, data_hz)
                        if x.size != filled_data.size and x.size > filled_data.size:
                            filled_data = np.append(filled_data, np.zeros(x.size - filled_data.size), axis=None)
                        elif x.size != filled_data.size and filled_data.size > x.size:
                            x = np.append(x, np.arange(q_T0_duration + data_hz,
                                                    q_T0_duration + data_hz * (filled_data.size - x.size), data_hz),
                                        axis=None)
                        if x.size != filled_data.size:
                            if x.size < filled_data.size:
                                filled_data = filled_data[0:x.size]
                            else:
                                x = x[0:filled_data.size]

                        f = interpolate.interp1d(x, filled_data, 'previous')
                        interpolated_data = f(q_T0)
                        signal_dict[signal_dict_key] = interpolated_data
                signal_dict["q_T0"] = q_T0

                # saving the converted signals to .mat file
                mat_filename = filename.split(sep=".")[0] + ".mat"
                scipy.io.savemat(path.join(temp_folder_path, mat_filename), signal_dict, appendmat=True,
                                long_field_names=True)
                extracted_files.append(mat_filename)
                logging.info("Done.")
                file_ID += 1
            except:
                print("Measurement file conversion failed for " + meas_file)

        rosbag_reader.terminate()

    def convert_mat(self):
        for mat_filepath in self.measurement_files:
            filename = mat_filepath.split(sep='\\')[-1]
            out_path = path.join(self.root, '_temp', filename)
            mat_file = scipy.io.loadmat(mat_filepath)
            mat_signal_list = list(mat_file.keys())
            mat_signal_list = mat_signal_list[3:]
            signal_list = [item[0] for item in self.config_signals]
            for mat_signal in mat_signal_list:
                # remove unnecessary signals from the mat file
                if mat_signal not in signal_list:
                    mat_file.pop(mat_signal)
                else:
                    for signals in self.config_signals:
                        if mat_signal in signals and len(signals) == 1:
                            continue
                        elif mat_signal in signals and len(signals) == 2:
                            mat_file[signals[-1]] = mat_file.pop(mat_signal)
            for config_signal in self.config_signals:
                # add missing signals to the mat file with value = 0
                if config_signal[0] not in mat_signal_list:
                    if len(config_signal) == 1:
                        mat_file[config_signal[0]] = np.zeros((1, 1))
                    else:
                        mat_file[config_signal[-1]] = np.zeros((1, 1))
            scipy.io.savemat(out_path, mat_file, appendmat=True, long_field_names=True)
