"""
Created on Wed Nov 28 12:18:46 2018
@author: jez2bp
"""
import os
import numpy
import logging
import scipy
import pandas
import numpy as np
import shutil
from _converter.conversion.reader.read_mdf import ReadMDF as mdf
from _converter.conversion.reader.read_mf4 import ReadMF4 as mf4
# from conversion.reader.read_d97 import ReadD97 as d97
from _converter.conversion.reader.read_d97v2 import ReadD97 as d97v2
from _converter.conversion.reader.read_plt import ReadPLT as plt
from _converter.conversion.expression import Expression
from _converter.conversion.dict2struct import CreateDictStruct
from scipy import interpolate
from scipy.io import savemat
from zipfile import ZipFile

class MainConversion():
    """The class implements functions and methods to drive the conversion from
    MDF file,to .mat file based on the .plt file
    """
    def __init__(self, settings, vehicle_profile, threadswitch, hookprogress=None):
        self.settings = settings
        self.vehicle_profile = vehicle_profile
        self.threadswitch = threadswitch
        self.hookprogress = hookprogress
        if self.settings["outputpath"] is None:
            self.settings["outputpath"] = ''
        self.outputlist = []
        self.evalnamelist = {}
        self.actualfile = None
        self.filesforread = self.proc_inputs(self.settings["inputpath"])

    def start(self):
        result = []
        if self.fileTypeOf(self.settings["pltpath"], ["plt"]):
            pltdict, sourcenames = self.get_plt_data(self.settings["pltpath"])
            pltcheck = self.check_plt_data(pltdict)
            if not self.is_invalid_pltdata(pltdict, sourcenames, pltcheck):
                fileindex = self.read_export_measurement(sourcenames, pltdict)
                self.setProgress(fileindex, set_to=100)
                if self.outputlist:
                    result = self.outputlist
                else:
                    logging.error("None of the files is convertable because of error or the current conversion options.")
            else:
                logging.error("Invalid .plt file")
        else:
            logging.error("This is not a .plt file!")
        return result

    def change_names_and_interpolate(self, signals, pltdict):
        """The function drives the signal renaming and the interpolation"""
        renamed_interpolated_dict = self.generate_new_signals_interp_and_rename(signals, pltdict)
        exphandler = Expression(pltdict, renamed_interpolated_dict)
        processed = exphandler.applicate_expression()
        return processed

    def rename_single_time_name(self, matfiles):
        """The method renames the time sample in the mat file"""
        for filename in matfiles:
            mat = scipy.io.loadmat(filename, mat_dtype=True)
            matsignalnames = mat.keys()
            timesignalname = 'timestamps'
            if 'time' in matsignalnames:
                timesignalname = 'time'
            mat['q_T0'] = mat[timesignalname]
            del mat[timesignalname]
            scipy.io.savemat(filename, mat, appendmat=True, long_field_names=True)

    def setProgress(self, fileindex, sectionnumber=0, set_to=None):
        if self.hookprogress is not None:
            if set_to is None:
                self.hookprogress(((fileindex+sectionnumber)/len(self.filesforread))*100)
            else:
                self.hookprogress(set_to)

    def remove_extension(self, abspath):
        filename = abspath
        if self.fileTypeOf(filename, ['dat', 'mf4', 'mdf', 'd97']):
            filenameparts = abspath.split('.')
            del filenameparts[-1]
            filename = '.'.join(filenameparts)
        return filename

    def get_data(self, sourcenames, pltdict, filename):
        """Drives the reading of the input files and the checking of the MDF
        signals in the extracted asmmdf.MDF() object
        """
        signals = False
        if self.fileTypeOf(filename, ['dat', 'mdf']):
            result = self.get_dat_mdf_data(filename, sourcenames, pltdict)
            if result['signals'] is not None and not result['skip']:
                signals = result['signals']
        elif self.fileTypeOf(filename, ['d97']):
            result = self.get_d97_data(filename, sourcenames, pltdict)
            if result['signals'] is not None and not result['skip']:
                signals = result['signals']
        elif self.fileTypeOf(filename, ['mf4']):
            result = self.get_mf4_data(filename, sourcenames, pltdict)
            if result['signals'] is not None and not result['skip']:
                signals = result['signals']
        else:
            logging.error("Unknown file extension: %s", filename)
        return signals

    def get_dat_mdf_data(self, filename, sourcenames, pltdict):
        """The function give back the needed asammdf.MDF() object and the list
        of the signals in the MDF() object in case of .dat and .mf4
        input files.
        """
        result = {'filename': filename, 'success': False,
                  'skip': False, 'signals': None,
                  'mdf_channels': None, 'error_msg': None}
        try:
            mdf_reader = mdf(filename, sourcenames)
            signals, mdf_channels = mdf_reader.read_channels()
            if mdf_channels:
                if signals is not None:
                    result['success'] = True
                    result['signals'] = signals
                    result['mdf_channels'] = mdf_channels
            else:
                result['success'] = False
                result['skip'] = True
                result['error_msg'] = "All signal is missing from the file"
        except Exception as readerror:
            logging.error("The reading was not succesfull! %s", str(readerror))
            result['success'] = False
            result['skip'] = True
            result['error_msg'] = readerror
        if result['success']:
            if self.settings["evals"] and self.settings["evalsignals"]:
                if not self.evaluate_mdf_signals(pltdict, mdf_channels):
                    result['success'] = False
                    result['skip'] = True
                    result['error_msg'] = 'File failed at signal evalution cheking'
            if not self.check_plt_signals_in_input(sourcenames, mdf_channels):
                if not self.settings["forcestate"]:
                    result['success'] = False
                    result['skip'] = True
                    result['error_msg'] = 'File failed because of missing PLT signal'
                else:
                    result['success'] = False
                    result['skip'] = False
                    result['error_msg'] = 'Missing PLT signal'
        return result

    def get_d97_data(self, filename, sourcenames, pltdict):
        """The function give back the needed asammdf.MDF() object and the list
        of the signals in the MDF() object in case of .d97 input files.
        """
        result = {'filename': filename, 'success': False,
                  'skip': False, 'mdf_object': None,
                  'mdf_channels': None, 'error_msg': None}
        try:
            d97parser = d97v2()
            d97parser.open(filename)
            signal_tuples = d97parser.getSignals(list(sourcenames))
            mdf_channels = []
            for signal in signal_tuples:
                mdf_channels.append(signal[0])
            if signal_tuples is not None:
                result['success'] = True
                result['signals'] = signal_tuples
                result['mdf_channels'] = mdf_channels

        except Exception as readerror:
            result['success'] = False
            result['skip'] = True
            result['error_msg'] = readerror
        if result['success']:
            if self.settings["evals"] and self.settings["evalsignals"]:
                if not self.evaluate_mdf_signals(pltdict, mdf_channels):
                    result['success'] = False
                    result['skip'] = True
                    result['error_msg'] = 'File failed at signal evalution cheking'
            if not self.check_plt_signals_in_input(sourcenames, mdf_channels):
                if not self.settings["forcestate"]:
                    result['success'] = False
                    result['skip'] = True
                    result['error_msg'] = 'File failed because of missing PLT signal'
                else:
                    result['success'] = False
                    result['skip'] = False
                    result['error_msg'] = 'Missing PLT signal'
        return result

    def get_mf4_data(self, filename, sourcenames, pltdict):
        """The function give back the needed signals Py list and the list
        of the signals in the Py list in case of .mf4 input files.
        """
        result = {'filename': filename, 'success': False,
                  'skip': False, 'signals': None,
                  'mdf_channels': None, 'error_msg': None}
        try:
            try:
                mf4_reader = mf4(filename, sourcenames)
                signals, mdf_channels = mf4_reader.read_channels()
            except:
                logging.warning('The new MF4READER failed therefore the measurement will be processed with the old ASAMMDF')
                logging.warning('Please send the PLT and MF4 files to the developers to investigate the issue and fix it!')
                mdf_reader = mdf(filename, sourcenames)
                signals, mdf_channels = mdf_reader.read_channels()
            if mdf_channels:
                if signals is not None:
                    result['success'] = True
                    result['signals'] = signals
                    result['mdf_channels'] = mdf_channels
            else:
                result['success'] = False
                result['skip'] = True
                result['error_msg'] = "All signal is missing from the file"
        except Exception as readerror:
            logging.error("The reading was not succesfull! %s", str(readerror))
            result['success'] = False
            result['skip'] = True
            result['error_msg'] = readerror
        if result['success']:
            if self.settings["evals"] and self.settings["evalsignals"]:
                if not self.evaluate_mdf_signals(pltdict, mdf_channels):
                    result['success'] = False
                    result['skip'] = True
                    result['error_msg'] = 'File failed at signal evalution cheking'
            if not self.check_plt_signals_in_input(sourcenames, mdf_channels):
                if not self.settings["forcestate"]:
                    result['success'] = False
                    result['skip'] = True
                    result['error_msg'] = 'File failed because of missing PLT signal'
                else:
                    result['success'] = False
                    result['skip'] = False
                    result['error_msg'] = 'Missing PLT signal'
        return result

    def check_plt_signals_in_input(self, sourcenames, mdf_channels):
        """The function returns True if the mdf_channels contains all of the
        of the plt signals else returns False and prints out the missing signals
        """
        if sorted(list(set(sourcenames))) != sorted(list(set(mdf_channels))):
            missingsignals = list(set(sourcenames)-set(mdf_channels))
            if missingsignals:
                if self.settings["forcestate"]:
                    logging.warning("Missing signals: %s", ",".join(missingsignals[:51]))
                else:
                    logging.error("Missing signals: %s", ",".join(missingsignals[:51]))
        return sorted(list(set(sourcenames))) == sorted(list(set(mdf_channels)))

    def channel_generator_without_master(self, mdf_object):
        channel_generator = mdf_object.iter_channels(skip_master=True)
        return channel_generator

    def proc_inputs(self, filename):
        result = []
        if os.path.isfile(filename):
            result.append(filename)
        else:
            foldername = os.path.dirname(filename)
            if os.path.isdir(foldername):
                result = sorted(self.listfiles(foldername, [".mf4", ".d97", ".dat", ".mdf", ".zip"]))
        return result

    def listfiles(self, path, extensions=[]):
        result = []
        with os.scandir(path) as files:
            for entry in files:
                parts = os.path.splitext(entry.name)
                if not extensions or parts[-1].lower() in list(map(lambda x: x.lower(), extensions)):
                    result.append(os.path.join(path, entry.name))
        return result

    def generate_new_signals_interp_and_rename(self, signals, pltdict):
        """The function implements the signal renaming and interpolation"""
        result = {}
        x_origin_list = []
        for timeitem in signals:
            x_origin_list.append(timeitem[2])
        minx, maxx = self.det_interpolation_limits(x_origin_list)
        for item in signals:
            x_values = item[2]
            y_values = item[1]
            try:
                x_g, y_g = self.cleanse_and_interpolate(x_values,
                                                                        y_values,
                                                                        minx,
                                                                        maxx,
                                                                        self.settings["samplingrate"]/1000)
            except ValueError as verror:
                logging.error("Interpolation is not possible with signal: %s, %s", str(item[0]), verror)
            sourcekey = []
            for target in pltdict.keys():
                sourcekey = pltdict[target]['name']
                if item[0] == sourcekey:
                    result[target] = y_g
        result['q_T0'] = x_g
        return result

    def det_interpolation_limits(self, x_origin_list):
        """The function determines the boundaries of the interpolation"""
        minx = None
        maxx = None
        if x_origin_list:
            mins = []
            maxx = []
            for x_origin in x_origin_list:
                try:
                    mins.append(min(x_origin))
                    maxx.append(max(x_origin))
                except ValueError as valerror:
                    logging.warning("Time limit determination failed %s", str(valerror))
            minx = max(mins)
            maxx = min(maxx)
        return minx, maxx

    def cleanse_and_interpolate(self, x_origin, y_origin, minx, maxx, delta_x):
        """
            The function erase the bad data from the signals and do the interpolation
        """
        data = {'x_origin': x_origin,
                'y_origin': y_origin}
        data_frame = pandas.DataFrame(data)
        data_frame_cleansed = data_frame.fillna(method='ffill')
        x_cleaned = data_frame_cleansed['x_origin']
        y_cleaned = data_frame_cleansed['y_origin']
        function = interpolate.interp1d(x_cleaned, y_cleaned, 'previous')
        x_new = np.arange(minx, maxx, delta_x)
        y_new = function(x_new)
        return x_new, y_new

    def export_mdf_to_mat(self, processed):
        """The method drives the export of the mdf object into mat"""
        # If the input path is a directory and the output path is a directory
        # or the output is not specfied we want to keep the original file name
        # of the MDF or D97 and save it to the output folder or the input
        # folder
        if (os.path.isdir(self.settings["inputpath"])
            and (os.path.isdir(self.settings["outputpath"])
                 or (self.settings["outputpath"] == ''))):
            self.export_folder_not_specified(processed)
        # If the input path is directory and the output path is not directory
        # but specified we want to name all of the input files as it is
        # specified. If we have more than one file in the input directory we
        # should indexing the names
        elif (os.path.isdir(self.settings["inputpath"])
              and (not os.path.isdir(self.settings["outputpath"])
                   and (self.settings["outputpath"] != ''))):
            self.export_folder_specified(processed)
        # If the input path is a file and the output path is a directory or
        # not specified we want to keep the original name of the file and save
        # it to the output folder or the input folder
        elif (os.path.isfile(self.settings["inputpath"])
              and (os.path.isdir(self.settings["outputpath"])
                   or (self.settings["outputpath"] == ''))):
            self.export_file_not_specified(processed)
        # If the input path is a file and the output path is not a directory
        # but specified we want to save the file to the folder of the output
        # with the specified name
        elif (os.path.isfile(self.settings["inputpath"])
              and (not os.path.isdir(self.settings["outputpath"])
                   and (self.settings["outputpath"] != ''))):
            self.export_file_specified(processed)

    def export_folder_not_specified(self, processed):
        filename = os.path.basename(self.actualfile['name'])
        if (os.path.isdir(self.settings["outputpath"])
                and (self.settings["outputpath"] != '')):
            dirname = self.settings["outputpath"]
        else:
            dirname = self.settings["inputpath"]
        exportname = os.path.join(dirname, filename)+'.mat'
        if "CAN_only" in self.vehicle_profile:
            processed = self.calc_missing_signals_CAN_only(processed)
        elif "DASy" in self.vehicle_profile:
            processed = self.calc_missing_signals_DASy(processed)
        elif "ADMA" in self.vehicle_profile:
            processed = self.calc_missing_signals_ADMA(processed)
        savemat(exportname, processed, long_field_names = True)
        self.outputlist.append(exportname)

    def export_folder_specified(self, processed):
        filename = os.path.basename(self.settings["outputpath"])
        if os.path.splitext(filename)[1].lower() != '.mat':
            if len(self.filesforread) > 1:
                filename = filename + str(self.actualfile['index']+1) + '.mat'
            else:
                filename = filename + '.mat'
        else:
            if len(self.filesforread) > 1:
                filename = os.path.splitext(filename)[0] + str(self.actualfile['index']+1) + os.path.splitext(filename)[1]
        dirname = os.path.dirname(self.settings["outputpath"])
        exportname = os.path.join(dirname, filename)
        if "CAN_only" in self.vehicle_profile:
            processed = self.calc_missing_signals_CAN_only(processed)
        elif "DASy" in self.vehicle_profile:
            processed = self.calc_missing_signals_DASy(processed)
        elif "ADMA" in self.vehicle_profile:
            processed = self.calc_missing_signals_ADMA(processed)
        savemat(exportname, processed, long_field_names = True)
        self.outputlist.append(exportname)

    def export_file_not_specified(self, processed):
        filename = os.path.basename(self.actualfile['name'])
        if (os.path.isdir(self.settings["outputpath"])
                and (self.settings["outputpath"] != '')):
            dirname = self.settings["outputpath"]
        else:
            dirname = os.path.dirname(self.settings["inputpath"])
        exportname = os.path.join(dirname, filename)+'.mat'
        if "CAN_only" in self.vehicle_profile:
            processed = self.calc_missing_signals_CAN_only(processed)
        elif "DASy" in self.vehicle_profile:
            processed = self.calc_missing_signals_DASy(processed)
        elif "ADMA" in self.vehicle_profile:
            processed = self.calc_missing_signals_ADMA(processed)
        savemat(exportname, processed, long_field_names = True)
        self.outputlist.append(exportname)

    def export_file_specified(self, processed):
        filename = os.path.basename(self.settings["outputpath"])
        if os.path.splitext(filename)[1].lower() != '.mat':
            filename = filename + '.mat'
        dirname = os.path.dirname(self.settings["outputpath"])
        exportname = os.path.join(dirname, filename)
        if "CAN_only" in self.vehicle_profile:
            processed = self.calc_missing_signals_CAN_only(processed)
        elif "DASy" in self.vehicle_profile:
            processed = self.calc_missing_signals_DASy(processed)
        elif "ADMA" in self.vehicle_profile:
            processed = self.calc_missing_signals_ADMA(processed)
        savemat(exportname, processed, long_field_names = True)
        self.outputlist.append(exportname)

    '''
    This method converts the left and right indicators of the jp data into the desired format
    In the JP format the 'Left_Index' and 'Right_Index' columns are the same with a 1 indicating the 'Left_Index' 
    and a 2 indicating a 'Right_Index'. In the desired format we have separate columns with a 1 indicating the corresponding index.
    '''
    def get_indicators(self, Left_Index):
        indicator_old = Left_Index
        Left_Index = []
        Right_Index = []
        for element in indicator_old:
            if element == 1:
                Left_Index.append(1)
                Right_Index.append(0)
            elif element == 2:
                Left_Index.append(0)
                Right_Index.append(1)
            else:
                Left_Index.append(0)
                Right_Index.append(0)
        return Left_Index, Right_Index

    def calc_missing_signals_DASy(self, processed):
        processed["AccelerationY_ESP"] = processed["AccelerationY_ESP"] / 2048
        processed["AccelerationX_ESP"] = processed["AccelerationX_ESP"] / 2048
        processed["yawRateESP"] = processed["yawRateESP"] / 16384
        processed["VelocityX_ESP"] = processed["VelocityX_ESP"] / 256
        processed["Left_Index"], processed["Right_Index"] = self.get_indicators(processed["Left_Index"])

        return processed

    def calc_missing_signals_ADMA(self, processed):
        processed["AccelerationY_ESP"] = processed["AccelerationY_ESP"] * 9.81
        processed["AccelerationX_ESP"] = processed["AccelerationX_ESP"] * 9.81
        processed["yawRateESP"] = processed["yawRateESP"] * 0.01745329

        return processed

    def calc_missing_signals_CAN_only(self, processed):
        # creating LUT for c1_dx, c2_dx
        x = np.array([0.000, 10.000, 50.000, 85.000, 100.000])
        y = np.array([0.001, 0.100, 0.400, 1.000, 1.000])
        LUT = interpolate.interp1d(x, y, kind='linear', fill_value='extrapolate')

        processed["c1_dx"] = LUT(processed["c1_dx"])
        processed["c2_dx"] = LUT(processed["c2_dx"])

        YawRate_sign = processed["yawRateESP_sign"]
        YawRate_sign = np.where(YawRate_sign == 0, YawRate_sign, -1)
        YawRate_sign = np.where(YawRate_sign < 0, YawRate_sign, 1)

        curve_norm_const = 0.00001525878906
        processed["c1"] = np.divide(np.multiply(processed["c1_orientation"], processed["c1_dx"]) +
                                    np.multiply(processed["c2_orientation"], processed["c2_dx"]), processed["c1_dx"] + processed["c2_dx"])
        processed["c2"] = np.divide(np.multiply(curve_norm_const * processed["c1_curvature"], processed["c1_dx"]) +
                                    np.multiply(curve_norm_const * processed["c2_curvature"], processed["c2_dx"]), processed["c1_dx"] + processed["c2_dx"])

        processed["SteeringAngle"] = processed["SteeringAngle"] * np.pi / 180
        processed["yawRateESP"] = np.multiply(processed["yawRateESP_deg_sec"] * np.pi / 180, YawRate_sign)

        processed.pop('c1_dx', None)
        processed.pop('c2_dx', None)
        processed.pop('c1_curvature', None)
        processed.pop('c2_curvature', None)
        processed.pop('c1_orientation', None)
        processed.pop('c2_orientation', None)
        processed.pop('yawRateESP_deg_sec', None)
        processed.pop('yawRateESP_sign', None)

        return processed

    def check_availability_plt(self, evallist, signallist):
        logging.info("PLT signal evaluation checking...")
        missing_signals = []
        missing_eval_data = {}
        errorswitch = False
        for eval_type in evallist:
            eval_type_signals = evallist[eval_type]
            for signalname in eval_type_signals:
                if signalname not in signallist:
                    errorswitch = True
                    missing_signals.append(signalname)
            if missing_signals:
                missing_eval_data[eval_type] = missing_signals
            missing_signals = []
        if errorswitch:
            for eval_mode in missing_eval_data:
                if missing_eval_data[eval_mode]:
                    logging.info("Missing signal from PLT to " + eval_mode + " evaluation: ")
                    logging.info(",".join(missing_eval_data[eval_mode]))
            return False
        logging.info('completed')
        return True

    def check_availability_mdf(self, evallist, signallist):
        logging.info("Input file signal checking...")
        missing_signals = []
        missing_eval_data = {}
        errorswitch = False
        for eval_type in evallist:
            eval_type_signals = evallist[eval_type]
            for signalname in eval_type_signals:
                if signalname not in signallist:
                    errorswitch = True
                    missing_signals.append(signalname)
            if missing_signals:
                missing_eval_data[eval_type] = missing_signals
            missing_signals = []
        if errorswitch:
            for eval_mode in missing_eval_data:
                if missing_eval_data[eval_mode]:
                    logging.info("Missing signal from PLT to " + eval_mode + " evaluation: ")
                    logging.info(",".join(missing_eval_data[eval_mode]))
            return False
        logging.info('completed')
        return True

    def evaluate_plt_signals(self, pltdict):
        """Prepares the evalution checking of the PLT signals"""
        result = {}
        pltnamelist = pltdict.keys()
        if self.settings["evalsignals"]:
            for eval_option_name, eval_option_data in self.settings["evalsignals"].items():
                self.evalnamelist[eval_option_name] = []
                for item in eval_option_data:
                    self.evalnamelist[eval_option_name].append(item['name'])
            result = self.check_availability_plt(self.evalnamelist,
                                                 list(set(pltnamelist)))
        return result

    def evaluate_mdf_signals(self, pltdict, mdf_channels):
        """Prepares the evalution checking of the MDF signals"""
        res = None
        if self.settings["evalsignals"]:
            mdfsignalnames = []
            for key in pltdict.keys():
                for channelname in mdf_channels:
                    if pltdict[key]['name'] == channelname:
                        mdfsignalnames.append(key)
            res = self.check_availability_mdf(self.evalnamelist,
                                              list(set(mdfsignalnames)))
        return res

    def fileTypeOf(self, filename, extensionlist):
        """The function returns True if the filename ends with the extension
           else returns False
        """
        return filename.split('.')[-1].lower() in list(map(lambda ext: ext.lower(), extensionlist))

    def is_invalid_pltdata(self, pltdict, sourcenames, check_result):
        """The function returns True if the result of the plt data is
           incorrect and plt data vs evaluation signal checking is
           incorrect"""
        if pltdict is None and sourcenames is None:
            return True
        if not check_result and self.settings["evalsignals"]:
            return True
        return False

    def get_plt_data(self, filename):
        """Drives the reading of the PLT file"""
        plt_handler = plt(filename)
        sourcenames, pltdict = plt_handler.extract_plt_data()
        return pltdict, sourcenames

    def check_plt_data(self, pltdict):
        """Drives the signal evaluation checking of the PLT data"""
        if self.settings["evals"]:
            return self.evaluate_plt_signals(pltdict)

    def extract_zipfile(self, fileindex):
        #in case zipfile is the input, extract it to temp folder
        os.mkdir("conv4KPI_tempdir")
        with ZipFile(self.filesforread[fileindex]) as zipObj:
            zipObj.extractall("conv4KPI_tempdir")
        files = sorted(self.listfiles("conv4KPI_tempDir", [".mf4", ".d97", ".dat", ".mdf"]))
        return files

    def read_measurement(self, files, fileindex, sourcenames, pltdict):
        filename = files[self.tempfileindex]
        self.setProgress(fileindex)
        logging.info("Reading %s ...", filename)
        self.actualfile = {'name': self.remove_extension(filename), 'index': self.subfileindex}
        signals = self.get_data(sourcenames, pltdict, filename)
        self.setProgress(fileindex, 1/3)
        return signals, filename

    def read_export_measurement(self, sourcenames, pltdict):
        # Reading the input files. If we got at least one processable file the
        # conversion will continues else the process will be stopped.
        fileindex = 0 #indexing the main files (raw meas and zipped files)
        self.subfileindex = 0 #indexing all the masurement files
        while not self.threadswitch.is_set() and fileindex < len(self.filesforread):
            files = []
            #check if zipped file, in that case unzip and return list of measurement files
            if os.path.splitext(self.filesforread[fileindex])[1].lower() == ".zip":
                files = self.extract_zipfile(fileindex)
            else:
                files.append(self.filesforread[fileindex])
            self.tempfileindex = 0 #indexing the files in zipped file
            while self.tempfileindex < len(files):
                #read measurement
                if not self.threadswitch.is_set():
                    signals, filename = self.read_measurement(files, fileindex, sourcenames, pltdict)
                    if not signals:
                        logging.warning("%s skipped", filename)
                        self.tempfileindex += 1
                        self.subfileindex += 1
                        continue
                # Renaming the signals based on PLT, interpolate signals, apply expressions.
                if not self.threadswitch.is_set():
                    logging.info('interpolation has been started')
                    processed = self.change_names_and_interpolate(signals, pltdict)
                    processed = self.rename_or_restructure_dict(processed)
                    logging.info('completed')
                    self.setProgress(fileindex, 2/3)
                    if not processed:
                        self.tempfileindex += 1
                        self.subfileindex += 1
                        continue
                # Export processed MDF object
                if not self.threadswitch.is_set():
                    logging.info('exporting has been started')
                    self.export_mdf_to_mat(processed)
                    logging.info('completed')
                    self.setProgress(fileindex, 1)
                self.tempfileindex += 1
                self.subfileindex += 1
            if os.path.splitext(self.filesforread[fileindex])[1].lower() == ".zip":
                shutil.rmtree("conv4KPI_tempDir")
            fileindex += 1
        return fileindex

    def rename_or_restructure_dict(self, processed):
        newDict = {}
        if not self.settings["openloop"]:
            for key in processed.keys():
                if (processed[key].dtype == numpy.float32 or processed[key].dtype == numpy.float16):
                    processed[key] = processed[key].astype(numpy.float64)
                temp_key = key
                if len(key) > 63:
                    key = key[:63]
                newDict[key.replace('[', '_').replace(']', '_').replace('.', '_')] = processed[temp_key]
        else:
            names = list(processed.keys())
            values = list(processed.values())
            for i in range(len(names)):
                names[i] = names[i].replace('[', '_').replace(']', '_')
                if (values[i].dtype == numpy.float32 or values[i].dtype == numpy.float16):
                    values[i] = values[i].astype(numpy.float64)
            dicthandle = CreateDictStruct(names, values)
            newDict = dicthandle.create_nested_dict()
        return(newDict)
