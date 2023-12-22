#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 11:56:15 2019
@author: jez2bp
"""

import threading
import logging
import libs.options
import libs.processes
import libs.config


class Application():

    def __init__(self, vehicle_profile):
        self.vehicle_profile = vehicle_profile
        options = libs.options.Options()
        settings = options.parseOptions()
        self.startConsole(settings)

    def startConsole(self, settings):
        logging.basicConfig(format='%(levelname)s: %(message)s', level=settings["loglevel"])
        config = libs.config.Config()
        threadswitch = threading.Event()
        threadswitch.clear()
        converter = libs.processes.ConvertProcesses(config, settings, self.vehicle_profile, threadswitch)
        converter.run()

