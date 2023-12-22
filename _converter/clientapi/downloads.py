#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 10:35:56 2018

@author: ote2bp
"""

import _converter.clientapi.sessions
import requests

class downloads(_converter.clientapi.sessions.Sessions):

    #get a download's version
    def checkversion(self,name):
        self.lastanswer = requests.post(self.apiurl + "/downloads/checkversion",
                                        data = {"name":name,
                                                "sessionid":self.sessionid})
        if self.lastanswer.status_code == 200:
            result = self.lastanswer.json()
            #we should use the resultcodes here, but have to move to the clientapi submodule
            if result["result"] == 0:
                return result["build"]
            else:
                #on error use the self.lastanswer te examine what happened
                return None
        else:
            #on error use the self.lastanswer te examine what happened
            return None
