#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:38:00 2019

@author: jez2bp
"""

import sys
import getopt
import logging
import os
from pathlib import Path
import json

class Options:
    """ The Option class handles the GUI mode flag or the running option
    arguments in CLI mode
    """

    def __init__(self):
        pass

    def parseOptions(self):
        file_sep = "/"
        plt_path = os.path.join(Path(os.path.dirname(os.path.abspath(__file__))).parent, 'tests', 'testPLT.PLT')
        plt_path = plt_path.split(sep="\\")
        plt_path = file_sep.join(plt_path)
        with open(os.path.join(Path(os.path.dirname(os.path.abspath(__file__))).parent.parent, '_temp',
                               'config.json')) as f:
            config = json.load(f)
        input_path = config["measurement_files"][0].split(sep='\\')
        input_path = file_sep.join(input_path[:-1]) + "/"
        output_path = os.path.join(Path(os.path.dirname(os.path.abspath(__file__))).parent.parent, '_temp')
        output_path = output_path.split(sep="\\")
        output_path = file_sep.join(output_path) + "/"
        result = {
            "pltpath": plt_path,
            "inputpath": input_path,
            "outputpath": output_path,
            "evals": [],
            "forcestate": True,
            "uploadstate": False,
            "loginstate": False,
            "samplingrate": 10,
            "loglevel": logging.INFO,
            "openloop": False
        }
        return result

    def help(self):
        helpstr = """\
The Conv4KPI convert .mf4 and .d97 files into .mat file based on the given .plt file.
If you want you can upload automatically the result files into the Online KPI Tool as input.
http://kpi.apps.de1.bosch-iot-cloud.com\n
Parameters:
    Mandatory parameters:
        --plt: the location of the plt file (filename with path)
        --input: the location of the MDF file (folder or file path). If a folder was given, all the .mf4 files will be the input what is in the folder
        --output: the location of the output .mat file (filename with path)

    Optional parameters:
        --force: force mode allows to the user if the MDF file not contains all signals what is on the PLT file or the signal is damaged, the conversion will continues with ignoring this data
        --upload: after the conversion the user can upload the generated .mat file into the Online KPI Tool. Till the authentication token expires not necessary to login into the Online KPI Tool account, the converted files will be uploaded into this account by default.
        --login: resets the upload profile settings. With this tag you need to login with username and password, if the authentication was successful, for the further uploads this account will be the default until the authentictaion token expires.
        --eval: checking the available signals based on the PLT and the measurement file for the evaluation in the Online KPI Tools
        --samplingrate: sampling rate of the result signal in ms.
            Options:
                aldhighspeed
                alc
                long
                sqt
                aldlowspeed
                mdi

Example commands:
conv4kpi.exe --plt=myfolder/myplt.plt --input=myfolder/ --output=result.mat

Basic command with all of the mandatory parameters. Based on the given .plt converts all .mf4 file into .mat what is in the myfolder/. If more than one file is given the output names will be result1.mat, result2.mat ... resultn.mat
conv4kpi.exe --plt=myfolder/myplt.plt --input=myfolder/mysourcefile.mf4 --output=resultfolder/result.mat
Converts the mysourcefile.mf4 to the result.mat file, based on the given .plt
"""
        print(helpstr)
