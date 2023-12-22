# -*- coding: utf-8 -*-
#******************************************************************************
#
#   Filename:    api_testclient.py
#   Author:      Paul, Dominik (POU2Fe)
#   Description: Test environment to call tests on the server side by HTTP
#                calls, example from command line or simple menu
#                used for manual tests
#
#******************************************************************************

#%% Notes
# 

#%% imports
import os
import configparser
# import requests
#from libs.resultcodes import ResultCodes as rc

# import from client API package
#import dataqueue as d
from dataqueue import Dataqueue as d
#import clientapi.dataqueue as d



#%% read the configuration
configfile = "setup.test.cfg"  # always loads the test config file
rootpath = os.path.dirname(os.path.abspath( __file__ ))
config = configparser.ConfigParser()
config.read(rootpath + "/defaults.cfg")
config.read(rootpath + "/" + configfile)
print(configfile + " config file loaded")

#%% main programm
if __name__ == "__main__":
    # the different urls
    apiurl = "http://kpi.apps.de1.bosch-iot-cloud.com/api"
    storagedir = os.path.join(rootpath, "downloads" )# the local savepath
    uploaddir = os.path.join(rootpath, "uploadfiles")
    
    # create the client API object
    apiDataqueue = d(config, apiPath)
    
    # login
    userid = apiDataqueue.login("test", "kpitest123")
    
    # get info about the whole queue
    queueInfoDict, preNumFiles, dirsize = apiDataqueue.listData()
    
    # remove all queue elements
    for n in list(queueInfoDict.keys()):
        status = apiDataqueue.removeData(n)
        if status:
            print("deleted queueid: \t"+str(n))
        else:
            print("something went wrong")
            
    # get the number of files
    queueInfoDict, preNumFiles, dirsize = apiDataqueue.listData()
    
    # upload the new files
    apiDataqueue.pushData(os.path.join(uploaddir, "test.mat"))
    apiDataqueue.pushData(os.path.join(uploaddir, "test.mf4"))
    apiDataqueue.pushData(os.path.join(uploaddir, "test.py"))
    apiDataqueue.pushData(os.path.join(uploaddir, "goodtest.zip"))
    apiDataqueue.pushData(os.path.join(uploaddir, "badtest.zip"))
    
    # get the number of files
    queueInfoDict, numberFiles, dirsize = apiDataqueue.listData()
    print(str(numberFiles-preNumFiles)+" files were uploaded")
    
    # download and pop Data
    file1 = apiDataqueue.downloadData(1, storagedir)
    file2 = apiDataqueue.popData(2, storagedir)
    print("first file saved at:\t"+file1)
    print("second file saved at:\t"+file2)
    
    
    # logut at the end
    if apiDataqueue.logout():
        print("logged out successfully")
    
        