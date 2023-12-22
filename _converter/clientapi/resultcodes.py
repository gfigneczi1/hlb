#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 18:11:52 2018

@author: ote2bp
"""

import enum

class ResultCodes(enum.IntEnum):
    #abort code, never returns
    ABORT = -1

    #general result codes, from 0
    OK = 0
    AUTHFAIL = 1
    MISSINGAUTHPARAMS = 2
    SESSIONERROR = 3
    SESSIONEXPIRED = 4
    UNKNOWNMETHOD = 5
    MISSINGPARAM = 6
    PERMISSIONDENIED = 7
    OLPREVENTED = 8

    #dataqueue related result from 50
    QUEUENOTFOUND = 51
    NOUPLOADFILE = 53
    FILENOTFOUND = 54
    FILETYPENOTALLOWED = 55
    QUEUEIDERROR = 56
    INSUFFICIENTSTORAGE = 57
    MISSINGPROJECTDESCRIPTION = 58
    MISSINGREPORTRESULT = 59
    MISSINGREPORTTEMPLATE = 60

    #script related result codes from 100
    SCRIPTTYPENOTALLOWED = 100
    SCRIPTDOESNTEXISTS = 101
    PARAMETERTYPEFAIL = 102
    UNKNOWNPARAMETER = 103
    UNKNOWNSCRIPTOPERATION = 104
    STARTFAIL = 105
    PROCESSINFOERROR = 106
    BADINPUT = 107
    NOAVAILABLESCRIPT = 108

    #user handling related result codes from 150
    USERNOTFOUND = 150
    EMAILNOTUNIQUE = 151
    BADEMAILFORMAT = 152
    NAMENOTUNIQUE = 153
    GROUPNOTFOUND = 154
    ALREADYMEMBER = 155
    NOTMEMBER = 156


    #download handling result codes from 200
    ENTRYNOTFOUND = 201
