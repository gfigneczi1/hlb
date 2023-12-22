#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 10:26:57 2020

@author: ote2bp
"""
import logging

class QueueLogger(logging.Handler):

    def __init__(self, queue):
        super().__init__()
        self.queue = queue

    def emit(self, item):
        self.queue.put(item)
