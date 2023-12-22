#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 20:48:41 2020

@author: ote2bp
"""
import os
import configparser

class Config:
    configpath = "appdata/config.cfg"
    
    def __init__(self):
        self.config = configparser.ConfigParser()
        self.config.read(self.configpath)
        self.structure()

    def structure(self):
        if not self.config.has_section("network"):
            self.config.add_section('network')
        if not self.config.has_section("paths"):
            self.config.add_section('paths')
        if "token" not in self.config["network"]:
            self.config.set('network', 'token', '')
        if "apiurl" not in self.config["network"]:
            self.config.set('network', "apiurl", "http://kpi.apps.de1.bosch-iot-cloud.com/api")
        if "input" not in self.config["paths"]:
            self.config.set('paths', 'input', '')
        if "output" not in self.config["paths"]:
            self.config.set('paths', 'output', '')
        if "plt" not in self.config["paths"]:
            self.config.set('paths', 'plt', '')
            
    def get(self, section, name):
        return self.config[section][name]
    
    def set(self, section, name, value):
        self.config[section][name] = value

    def save(self):
        if not os.path.exists(os.path.dirname(self.configpath)):
            os.makedirs(os.path.dirname(self.configpath))
        with open(self.configpath, 'w') as configfile:
            self.config.write(configfile)
