# -*- coding: utf-8 -*-
# ******************************************************************************
#
#   Filename:    dataqueue.py
#   Author:      Paul, Dominik (POU2Fe)
#   Description: The client side dataqueue class
#
# ******************************************************************************

# %% Notes


# %% imports
import requests
from _converter.clientapi.sessions import Sessions
import os


# %% the client side dataqueue class
class Dataqueue(Sessions):
    '''
     The client side dataqueue class
     The following methods are included:
        queueid = self.pushData(filename):
            - uploads the given file
            - returns the new queueid
        savepath = self.downloadData(queueid, storageDir):
            - can request the download of the queue element with the queueid
              and store it at storageDir
            - returns the path where the file was saved
        queueids = self.popData(userdir, storageDir):
            - can request the download and removalof the queue element with
              the queueid and store it at storageDir
            - returns the path where the file was saved
        status = self.removeData(queueid):
            - removes the queue element defined by the userid and the queueid
            - returns whether it was removed successful or not with the boolean
              values True/False
        numFiles, totalSize = self.queueInfo():
            - return the total number of files and the total Size of the whole
              queue
        queueInfoDict, numberFiles, dirsize = self.listData():
            - returns the whole information about the queue of the user defined
              by the userid
            - queueInfoDict is a dict. Each element has a queueid as the key
              and the value is an dict as well.
              This dict contains until now the following information:
              {"last modified":modTime}
            - numberFiles is the number of files in the user's queue and
              dirsize is the size of the user's queue in bytes
    '''
    # %% Methods
    # push one element to the queue and return the queueid assigned to it
    # returns the status
    def pushData(self, filenames):
        '''
        queueid = self.pushData(filename):
            - uploads the given file
            - returns the new queueid
        '''
        if type(filenames).__name__ == 'string':
            filenames = [filenames]
        files = {}
        num = 1
        for filename in filenames:
            if os.path.isfile(filename):
                fp = open(filename, "rb")
                files["file"+str(num)] = fp
                num = num + 1
        self.lastanswer = requests.post(self.apiurl+"/dataqueue/push",
                                        data={"sessionid": self.sessionid},
                                        files=files)
        for num in files.keys():
            files[num].close()
        if self.lastanswer.status_code == 200:
            result = self.lastanswer.json()
            if result["result"] == 0:
                return result["queueid"]
            else:
                return None
        else:
            return None

    # download data as a zip file from the queue by its queueid
    # important: client has to check if the response contains an error in json
    #            format otherwise it contains the file/zipfile

    def downloadData(self, queueid, storageDir):
        '''
        savepath = self.downloadData(queueid, storageDir):
            - can request the download of the queue element with the queueid
              and store it at storageDir returns the path where the file was
              saved
        '''
        self.lastanswer = requests.post(self.apiurl+"/dataqueue/download",
                                        data={"sessionid": self.sessionid,
                                              "queueid": queueid})
        if self.lastanswer.status_code == 200:
            result = self.lastanswer
            # check if the downloaded data is a file or json -> json when file
            # not found
            # structure of ["Content-Type"] when its json: "applicaton/json"
            filetype = result.headers["Content-Type"].rsplit("/")[-1].lower()
            if filetype == "json":
                return None
            # if the response content is not json, a file is sended back
            else:
                # get the filename out of the Content-Disposition
                # structure of ["Content-Disposition"]:
                # "attachment:filename=<filename>
                filedescript = result.headers["Content-Disposition"]
                filename = filedescript.rsplit("=")[-1]
                # create the savepath
                savepath = os.path.join(storageDir, filename)
                # save the downloaded file
                with open(savepath, "wb") as f:
                    f.write(result.content)
                return (savepath)
        else:
            return None

    # download and remove an element from the queue by id
    def popData(self, queueid, storageDir):
        '''
        queueids = self.popData(userdir, storageDir):
            - can request the download and removalof the queue element with the
              queueid and store it at storageDir
            - returns the path where the file was saved
        '''
        self.lastanswer = requests.post(self.apiurl+"/dataqueue/pop",
                                        data={"sessionid": self.sessionid,
                                              "queueid": queueid})

        if self.lastanswer.status_code == 200:
            result = self.lastanswer
            # check if the downloaded data is a file or json -> json when file
            # not found
            # structure of ["Content-Type"] when its json: "applicaton/json"
            filetype = result.headers["Content-Type"].rsplit("/")[-1].lower()
            if filetype == "json":
                return None
            # if the response content is not json, a file is sended back
            else:
                # get the filename out of the Content-Disposition
                # structure of ["Content-Disposition"]:
                # "attachment:filename=<filename>
                filedescript = result.headers["Content-Disposition"]
                filename = filedescript.rsplit("=")[-1]
                # create the savepath
                savepath = os.path.join(storageDir, filename)
                # save the downloaded file
                with open(savepath, "wb") as f:
                    f.write(result.content)
                return (savepath)
        else:
            return None

    # remove an element from the queue by ID
    # returns the status
    def removeData(self, queueid):
        '''
        status = self.removeData(queueid):
            - removes the queue element defined by the userid and the queueid
            - returns whether it was removed successful or not with the boolean
              values True/False
        '''
        self.lastanswer = requests.post(self.apiurl+"/dataqueue/remove",
                                        data={"sessionid": self.sessionid,
                                              "queueid": queueid})
        if self.lastanswer.status_code == 200:
            result = self.lastanswer.json()
            if result["result"] == 0:
                return True
            else:
                return None
        else:
            return None

    # get information about the whole queue
    # returns the number of elements and file size summary
    def queueInfo(self):
        '''
        numFiles, totalSize = self.queueInfo():
            - return the total number of files and the total Size of the whole
              queue
        '''
        self.lastanswer = requests.post(self.apiurl+"/dataqueue/queueInfo",
                                        data={"sessionid": self.sessionid})
        if self.lastanswer.status_code == 200:
            result = self.lastanswer.json()
            if result["result"] == 0:
                # returns the answer of the server without the result element
                return (result["number of files"],
                        result["total size in bytes"])
            else:
                return None
        else:
            return None

    # get a list of the user's stored data elements
    # returns a array of all queue IDs, sum filesize, number of files and
    # creation datetime
    def listData(self):
        '''
        queueInfoDict, numberFiles, dirsize = self.listData():
            - returns the whole information about the queue of the user defined
              by the userid
            - queueInfoDict is a dict. Each element has a queueid as the key
              and the value is an dict as well.
              This dict contains until now the following information:
              {"last modified":modTime}
            - numberFiles is the number of files in the user's queue and
              dirsize is the size of the user's queue in bytes
        '''
        self.lastanswer = requests.post(self.apiurl+"/dataqueue/list",
                                        data={"sessionid": self.sessionid})
        if self.lastanswer.status_code == 200:
            result = self.lastanswer.json()
            if result["result"] == 0:
                # returns the answer of the server without the result element
                return (result["queueids info"],
                        result["number of files"],
                        result["total size in bytes"])
            else:
                return None
        else:
            return None
