#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 10:35:56 2018

@author: ote2bp
"""

import _converter.clientapi.sessions
import requests

class Processes(_converter.clientapi.sessions.Sessions):

    #query signals for a processline
    def listProcesssignals(self,plname):
        result = None
        self.lastanswer = requests.post(self.apiurl + "/process/processsignals",data = {"sessionid":self.sessionid,"plname":plname})
        if self.lastanswer.status_code == 200:
            answer = self.lastanswer.json()
            if answer["result"] == 0:
                result = answer["signals"]
        return result
