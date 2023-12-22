"""
Created on Wed Nov 28 12:18:46 2018
@author: jez2bp
"""
import os
import numpy
import asammdf
import logging
import scipy 
import pandas
import numpy as np
from conversion.reader.read_mdf import ReadMDF as mdf
from conversion.reader.read_d97 import ReadD97 as d97
from conversion.reader.read_plt import ReadPLT as plt
from conversion.expression import Expression
from scipy import interpolate

class MainConversion():
    """The class implements functions and methods to drive the conversion from
    MDF file,to .mat file based on the .plt file
    """
    def __init__(self, settings, threadswitch, hookprogress=None):
        self.settings = settings
        self.threadswitch = threadswitch
        self.hookprogress = hookprogress
        if self.settings["outputpath"] is None:
            self.settings["outputpath"] = ''
        self.outputlist = []
        self.evalnamelist = {}
        self.actualfile = None
        self.filesforread = self.proc_inputs(self.settings["inputpath"])

    def start(self):
        if self.fileTypeOf(self.settings["pltpath"], ["plt"]):
            pltdict, sourcenames = self.get_plt_data(self.settings["pltpath"])
            pltcheck = self.check_plt_data(pltdict)
            if not self.is_invalid_pltdata(pltdict, sourcenames, pltcheck):
                # Reading the input files. If we got at least one processable file the
                # conversion will continues else the process will be stopped.
                fileindex = 0
                while not self.threadswitch.is_set() and fileindex < len(self.filesforread):
                    #read measurement
                    if not self.threadswitch.is_set():
                        filename = self.filesforread[fileindex]
                        self.setProgress(fileindex)
                        logging.info("Reading %s ...", filename)
                        self.actualfile = {'name': self.remove_extension(filename), 'index': fileindex}
                        mdf_object = self.get_data(sourcenames, pltdict, filename)
                        self.setProgress(fileindex, 1/3)
                        if not mdf_object:
                            logging.warning("%s skipped", filename)
                            fileindex += 1
                            continue
                    # Renaming the signals based on PLT, interpolate signals, apply expressions.
                    if not self.threadswitch.is_set():
                        logging.info('interpolation has been started')
                        processed = self.change_names_and_interpolate(mdf_object, pltdict)
                        logging.info('completed')
                        self.setProgress(fileindex, 2/3)
                        if not processed:
                            fileindex += 1
                            continue
                    # Export processed MDF object
                    if not self.threadswitch.is_set():
                        logging.info('exporting has been started')
                        self.export_mdf_to_mat(processed)
                        logging.info('completed')
                        self.setProgress(fileindex, 1)
                    fileindex += 1
                self.setProgress(fileindex, set_to=100)
                if self.outputlist:
                    self.rename_single_time_name(self.outputlist)
                else:
                    logging.error("None of the files is convertable because of error or the current conversion options.")
            else:
                logging.error("Invalid .plt file")
        else:
            logging.error("This is not a .plt file!")

    def change_names_and_interpolate(self, shrinked_mdf, pltdict):
        """The function drives the signal renaming and the interpolation"""
        renamed_interpolated_mdf = self.generate_new_signals_interp_and_rename(shrinked_mdf, pltdict)
        exphandler = Expression(pltdict, renamed_interpolated_mdf)
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
        shrinked_mdf = False
        if self.fileTypeOf(filename, ['dat', 'mf4', 'mdf']):
            result = self.get_dat_mdf_data(filename, sourcenames, pltdict)
            if isinstance(result['mdf_object'], asammdf.MDF) and not result['skip']:
                shrinked_mdf = result['mdf_object']
        elif self.fileTypeOf(filename, ['d97']):
            result = self.get_d97_data(filename, sourcenames, pltdict)
            if isinstance(result['mdf_object'], asammdf.MDF) and not result['skip']:
                shrinked_mdf = result['mdf_object']
        else:
            logging.error("Unknown file extension: %s", filename)
        return shrinked_mdf

    def get_dat_mdf_data(self, filename, sourcenames, pltdict):
        """The function give back the needed asammdf.MDF() object and the list
        of the signals in the MDF() object in case of .dat and .mf4
        input files.
        """
        result = {'filename': filename, 'success': False,
                  'skip': False, 'mdf_object': None,
                  'mdf_channels': None, 'error_msg': None}
        try:
            mdf_reader = mdf(filename, sourcenames)
            mdfobject, mdf_channels = mdf_reader.read_channels()
            if mdf_channels:
                if isinstance(mdfobject, asammdf.MDF):
                    result['success'] = True
                    result['mdf_object'] = mdfobject
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
            d97parser = d97(filename)
            d97parser.D97Read()
            mdfobject = d97parser.get_signal(list(set(sourcenames)))
            mdf_channels = []
            for signal in self.channel_generator_without_master(mdfobject):
                mdf_channels.append(signal.name)
            if isinstance(mdfobject, asammdf.MDF):
                result['success'] = True
                result['mdf_object'] = mdfobject
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
                result = sorted(self.listfiles(foldername, [".mf4", ".d97", ".dat", ".mdf"]))
        return result

    def listfiles(self, path, extensions=[]):
        result = []
        with os.scandir(path) as files:
            for entry in files:
                parts = os.path.splitext(entry.name)
                if not extensions or parts[-1].lower() in list(map(lambda x: x.lower(), extensions)):
                    result.append(os.path.join(path, entry.name))
        return result

    def generate_new_signals_interp_and_rename(self, mdf_data, pltdict):
        """The function implements the signal renaming and interpolation"""
        result = asammdf.MDF()
        x_origin_list = []
        for timeitem in self.channel_generator_without_master(mdf_data):
            x_origin_list.append(timeitem.timestamps)
        minx, maxx = self.det_interpolation_limits(x_origin_list)
        for item in self.channel_generator_without_master(mdf_data):
            x_values = item.timestamps
            y_values = item.samples
            try:
                x_g, y_g = self.cleanse_and_interpolate(x_values,
                                                                        y_values,
                                                                        minx,
                                                                        maxx,
                                                                        self.settings["samplingrate"]/1000)
            except ValueError as verror:
                logging.error("Interpolation is not possible with signal: %s, %s", str(item.name), verror)
            sourcekey = []
            for target in pltdict.keys():
                sourcekey = pltdict[target]['name']
                if item.name == sourcekey:
                    gen_signal = [asammdf.mdf.Signal(samples=y_g, name=target, timestamps=x_g)]
                    result.append(gen_signal)
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
        processed = self.finalize_signaltypes(processed)
        filename = os.path.basename(self.actualfile['name'])
        if (os.path.isdir(self.settings["outputpath"])
                and (self.settings["outputpath"] != '')):
            dirname = self.settings["outputpath"]
        else:
            dirname = self.settings["inputpath"]
        exportname = os.path.join(dirname, filename)+'.mat'
        processed.export(fmt="mat", filename=exportname,
                         single_time_base=True)
        self.outputlist.append(exportname)

    def export_folder_specified(self, processed):
        processed = self.finalize_signaltypes(processed)
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
        processed.export(fmt="mat", filename=exportname,
                         single_time_base=True)
        self.outputlist.append(exportname)

    def export_file_not_specified(self, processed):
        processed = self.finalize_signaltypes(processed)
        filename = os.path.basename(self.actualfile['name'])
        if (os.path.isdir(self.settings["outputpath"])
                and (self.settings["outputpath"] != '')):
            dirname = self.settings["outputpath"]
        else:
            dirname = os.path.dirname(self.settings["inputpath"])
        exportname = os.path.join(dirname, filename)+'.mat'
        processed.export(fmt="mat", filename=exportname,
                         single_time_base=True)
        self.outputlist.append(exportname)

    def export_file_specified(self, processed):
        processed = self.finalize_signaltypes(processed)
        filename = os.path.basename(self.settings["outputpath"])
        if os.path.splitext(filename)[1].lower() != '.mat':
            filename = filename + '.mat'
        dirname = os.path.dirname(self.settings["outputpath"])
        exportname = os.path.join(dirname, filename)
        processed.export(fmt="mat", filename=exportname,
                         single_time_base=True)
        self.outputlist.append(exportname)

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

    def finalize_signaltypes(self, processed):
        typecasted = asammdf.MDF()
        for processed_signal in processed.iter_channels(skip_master=True):
            if (processed_signal.samples.dtype == numpy.float32 or processed_signal.samples.dtype == numpy.float16):
                processed_signal.samples = processed_signal.samples.astype(numpy.float64)
                typecasted.append([processed_signal])
            else:
                typecasted.append([processed_signal])
        return typecasted
