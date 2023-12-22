#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 11:56:01 2019
@author: jez2bp
"""

class Buildinfo:
    """
    The class implements functions to read and modify the buildnumber
    """
    def __init__(self):
        self.buildnumber = 115

    def read(self):
        return self.buildnumber
