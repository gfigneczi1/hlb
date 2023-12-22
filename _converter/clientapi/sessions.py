#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 10:36:30 2018

@author: ote2bp
"""
import requests

class Sessions:
    def __init__(self,config,apiurl):
        self.config = config
        self.apiurl = apiurl
        self.sessionid = None
        self.lastanswer = None

    #login to the service
    def login(self,username,password):
        result = None
        self.lastanswer = requests.post(self.apiurl + "/users/login",data = {"username":username, "password":password})
        if self.lastanswer.status_code == 200:
            answer = self.lastanswer.json()
            if answer["result"] == 0:
                self.sessionid = answer["sessionid"]
                result = answer["userid"]
        return result

    #login to the service with stored token
    def tokenLogin(self,token):
        result = None
        self.lastanswer = requests.post(self.apiurl + "/users/tokenlogin",data = {"token":token})
        if self.lastanswer.status_code == 200:
            answer = self.lastanswer.json()
            if answer["result"] == 0:
                self.sessionid = answer["sessionid"]
                result = answer["userid"]
        return result

    #when the user logged in, able to require a long session token to store this token instead of credentials
    #use the tokenLogin method to login with this token
    def requireToken(self):
        result = None
        self.lastanswer = requests.post(self.apiurl + "/users/requiretoken",data = {"sessionid":self.sessionid})
        if self.lastanswer.status_code == 200:
            answer = self.lastanswer.json()
            if answer["result"] == 0:
                result = answer["token"]
        return result

    #get stored information about logged user
    def userinfo(self):
        result = None
        self.lastanswer = requests.post(self.apiurl + "/users/info",data = {"sessionid":self.sessionid})
        if self.lastanswer.status_code == 200:
            answer = self.lastanswer.json()
            if answer["result"] == 0:
                result = answer["user"]
        return result

    #logout from online KPI
    def logout(self):
        result = None
        self.lastanswer = requests.post(self.apiurl + "/users/logout",data = {"sessionid":self.sessionid})
        if self.lastanswer.status_code == 200:
            result = self.lastanswer.json()
            if result["result"] == 0:
                self.sessionid = None
                result = True
        return result
