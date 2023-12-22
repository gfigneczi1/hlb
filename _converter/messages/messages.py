#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:58:54 2019

@author: jez2bp
"""

import logging
from _converter.appdata.buildinfo import Buildinfo

class Messages:

    def __init__(self, message_id, parameters=None):
        self.message_id = message_id
        self.parameters = parameters
        try:
            self.execute()
        except ValueError:
            pass

    def execute(self):
        if self.message_id == 'stop':
            self.print_stop()
        elif self.message_id == 'read':
            self.print_read()
        elif self.message_id == 'outputlist':
            self.print_outputlist()
        elif self.message_id == 'skip_message':
            self.print_skipmessage()
        elif self.message_id == 'fail_message':
            self.print_failmessage()
        elif self.message_id == 'unknown_file':
            self.print_unknown_extension()
        elif self.message_id == 'missing_signal':
            self.print_missing_signals()
        elif self.message_id == 'completed':
            self.print_completed()
        elif self.message_id == 'processing':
            self.print_processing()
        elif self.message_id == 'export':
            self.print_export()
        elif self.message_id == 'plt_signalcheck':
            self.print_plt_signalcheck()
        elif self.message_id == 'missing_plt_signal':
            self.print_missing_plt_signals()
        elif self.message_id == 'input_signalcheck':
            self.print_input_signalcheck()
        elif self.message_id == 'upload':
            self.print_upload_message()
        elif self.message_id == 'connection_required':
            self.print_connection_message()
        elif self.message_id == 'signal_eval_check':
            self.print_signal_eval_check()
        elif self.message_id == "authfailed":
            self.print_authfailed()

    def print_read(self):
        filepath = self.parameters
        logging.info("Reading %s ...", filepath)

    def print_processing(self):
        logging.info("Processing...")

    def print_export(self):
        logging.info("Exporting...")

    def print_stop(self):
        logging.info("Stopped")

    def print_completed(self):
        logging.info("Completed")

    def print_outputlist(self):
        createdfiles = self.parameters
        if createdfiles:
            for filename in createdfiles:
                logging.info("%s has been created", filename)
        else:
            logging.info("None of the files is convertable because of error or the current conversion options.")

    def print_skipmessage(self):
        logging.info("The file is not processable because of error / the current conversion option.")

    def print_failmessage(self):
        logging.info("Processing failed. Conversion is not possible.")

    def print_unknown_extension(self):
        filename = self.parameters
        logging.info("Uknown file extension:")
        logging.info(filename)

    def print_missing_signals(self):
        sourcenames, mdf_channels, force = self.parameters
        missingsignals = list(set(sourcenames)-set(mdf_channels))
        if missingsignals:
            if force:
                logging.warning("Missing signals: %s", ",".join(missingsignals[:51]))
            else:
                logging.error("Missing signals: %s", ",".join(missingsignals[:51]))

    def print_plt_signalcheck(self):
        logging.info("PLT SIGNAL EVALUATION CHECKING...")

    def print_missing_plt_signals(self):
        missing_signal_data = self.parameters
        for eval_mode in missing_signal_data:
            missing_signals = missing_signal_data[eval_mode]
            logging.info("Missing signal from PLT to " + eval_mode + " evaluation: ")
            missing_signal_text = ''
            for name in missing_signals:
                if name is not missing_signals[-1]:
                    missing_signal_text = missing_signal_text + name + ", "
                else:
                    missing_signal_text = missing_signal_text + name
            logging.info(missing_signal_text)

    def print_input_signalcheck(self):
        logging.info("INPUT FILE SIGNAL EVALUATION CHECKING...")

    def print_upload_message(self):
        logging.info("UPLOADING...")

    def print_connection_message(self):
        logging.info("You are using at least one option what requires cloud connection!")
        logging.info("Please login into the Online KPI Tool account or register...")
        logging.info("http://kpi.apps.de1.bosch-iot-cloud.com\n")

    def print_signal_eval_check(self):
        logging.info("REQUESTING SIGNAL NAMES FOR SIGNAL EVALUATION CHECKING...")

    def print_authfailed(self):
        optionname = self.parameters
        message = str(optionname) + " is not possible because authetication failed!"
        logging.info(message)

class ErrorMessages:

    def __init__(self, message_id, parameters=None):
        self.message_id = message_id
        self.parameters = parameters
        build = Buildinfo()
        self.build = build.read()
        self.execute()

    def execute(self):
        if self.message_id == 'upload':
            self.error_upload()
        elif self.message_id == 'token_auth_failed':
            self.error_token_auth()
        elif self.message_id == 'build_outdated':
            self.error_build_outdated()
        elif self.message_id == 'missing_password':
            self.error_missing_password()
        elif self.message_id == 'failed_auth_password':
            self.error_failed_auth_password()
        elif self.message_id == 'connection_not_possile':
            self.error_connection()
        elif self.message_id == 'invalid_sampling_rate':
            self.error_invalid_sampling_rate()
        elif self.message_id == 'invalid_sampling_type':
            self.error_invalid_sampling_type()
        elif self.message_id == 'read_error':
            self.error_reading()
        elif self.message_id == 'asammdf_readerror':
            self.error_asammdf_reading()
        elif self.message_id == 'error_memory_limit':
            self.error_reached_memory_limit()
        elif self.message_id == 'corrupted_signal':
            self.error_corrupted_signal()
        elif self.message_id == 'expression_error':
            self.error_expression()
        elif self.message_id == "limit_error":
            self.error_limit()
        elif self.message_id == 'nodir':
            self.error_nodir()
        elif self.message_id == 'invalid_samplingrate':
            self.error_wrong_samplingrate()

    def error_upload(self):
        current_dataqueue = self.parameters
        logging.warning("Something went wrong while uploading")
        logging.warning("Status code: " + str(current_dataqueue.lastanswer.status_code))
        logging.warning('Maybe the file is too big!')

    def error_token_auth(self):
        logging.error("AUTHENTICATION FAILED WITH TOKEN")

    def error_build_outdated(self):
        logging.warning("Your current Conv4KPI build is outdated! Please  use the newest version!\nDownload: http://kpi.apps.de1.bosch-iot-cloud.com")

    def error_missing_password(self):
        logging.warning("Password is missing!")

    def error_failed_auth_password(self):
        logging.error("AUTHENTICATION FAILED! Try again!")

    def error_connection(self):
        logging.error('CONNECTION IS NOT POSSIBLE!')

    def error_interpolation(self):
        name, error = self.parameters
        logging.warning("Interpolation is not possible with signal: "+str(name))
        logging.error("Problem: " + str(error))

    def error_invalid_sampling_rate(self):
        logging.error("Invalid sampling rate!: The value need to be between 1 and 1000")

    def error_invalid_sampling_type(self):
        logging.error("Invalid sampling type! Sampling rate need to be an integer!")

    def error_reading(self):
        error = self.parameters
        message = "The reading was not succefull! " + str(error)
        logging.error(message)

    def error_asammdf_reading(self):
        error = self.parameters
        message = "Error while reading the input file! " + str(error)
        logging.warning(message)

    def error_reached_memory_limit(self):
        logging.error("Memory error!")

    def error_corrupted_signal(self):
        channel, groupindex, channelindex, signalerror = self.parameters
        logging.error("Signal is corrupted:")
        logging.warning("Signal name:" + str(channel))
        logging.warning("Group index:" + str(groupindex))
        logging.warning("Channel index:" + str(channelindex))
        logging.error(str(signalerror))

    def error_option(self):
        logging.error("Invalid option! %s", str(self.parameters))

    def error_expression(self):
        signalname, expression = self.parameters
        logging.warning("UNKNOWN EXPRESSION")
        logging.warning('SIGNAL NAME: ' + str(signalname))
        logging.warning('EXPRESSION: ' + str(expression))

    def error_limit(self):
        valerror = self.parameters
        message = "Time limit determination failed" + str(valerror)
        logging.warning(message)

    def error_nodir(self):
        directory = self.parameters
        message = str(directory) + " not exists"
        logging.error(message)

    def error_wrong_samplingrate(self):
        value, minvalue, maxvalue = self.parameters
        message_1 = "The value of the sampling rate ({value}) must be a number between {minvalue} and {maxvalue} ms.".format(value=value, minvalue=minvalue, maxvalue=maxvalue)
        message_2 = "The program running with the default 5 ms value."
        logging.warning(message_1)
        logging.warning(message_2)
