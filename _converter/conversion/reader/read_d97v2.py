#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 15:02:30 2021

@author: ote2bp
"""
import re
import struct
import numpy as np

"""
    Bosch D97 file reader class
    Missing features:
        - handling signal's TPOS field
        - handling bit fields
        - handling ASCII float format
        - recognize file as a valid D97
"""

class ReadD97:
    # For section position extracting
    SECTIONEXP = re.compile(r'(\[.*\])')
    CONTROLEXP = re.compile(r'\[CONTROL\]', re.IGNORECASE)
    SIGNALEXP = re.compile(r'\[\SIGNAL(\d+)\]', re.IGNORECASE)
    DATAEXP = re.compile(r'\[DATA\]', re.IGNORECASE)

    """
    ********* Public methods start here ********
    """
    #opens a D97 file
    def open(self, filename):
        self.sectionmap = {
            'CONTROL': {},
            'SIGNALS': {},
            'SIGNALINDEX' : {},
            #signal or framemap, depends on layout(DATAMATRIX)
            'MAP' : {},
            'TIME': {}
            }
        self.dataposition = 0
        self.fp = open(filename, 'rb', buffering=1024*1024)
        self._buildSectionmap()
        self._buildMap()

    #closes a D97 file
    def close(self):
        self.fp.close()

    #filters signalnames from d97 using regexp expressions
    def signalnames(self, refilter):
        if refilter == "*":
            refilter = ".+"
        rexp = re.compile(refilter)
        return list(filter(lambda signalname: rexp.fullmatch(signalname) is not None, self.sectionmap["SIGNALINDEX"].values()))

    #gets signal data by signalname list, the result is a list of tuples [(name, samples, timesamples)]
    def getSignals(self, signallist):
        result = []
        for signalname in signallist:
            signal = self.getSignal(signalname)
            if signal is not None:
                result.append(signal)
        return result

    #gets a signal data by name, the result is a tuple (name, samples, timesamples)
    def getSignal(self, signalname):
        result = None
        if signalname in self.sectionmap["SIGNALS"].keys():
            signalinfo = self.sectionmap["SIGNALS"][signalname]
            #reads the signals frame by frame, every frame contains a repetition length of signaldata, on DPOS position
            if self.sectionmap["CONTROL"]["DATAMATRIX"] == "FRAME_BY_*":
                result = self._getSignal_FrameBy(signalname, signalinfo)
            #reads the signals sequential
            elif self.sectionmap["CONTROL"]["DATAMATRIX"] == "SV_BY_*":
                result = self._getSignal_SvBy(signalname, signalinfo)
        return result
    """
    ********* Public methods end here ********
    ********* Low level methods start here *******
    """
    #reads only one signal, the result is a tuple (name, samples, timesamples)
    def _getSignal_FrameBy(self, signalname, signalinfo):
        #first dataposition in a frame
        dataposition = int(signalinfo["DPOS"])*8
        samplecount = int(signalinfo["NEXTDPOS"]) - int(signalinfo["DPOS"])
        #number of samples in a frame
        repetition = int(signalinfo["RPT"])
        samples = self._readFrameBy(dataposition, repetition, samplecount)
        timesamples = self._getTime(signalinfo)
        return (signalname, samples, timesamples)

    #reads only one signal, the result is a tuple (name, samples, timesamples)
    def _getSignal_SvBy(self, signalname, signalinfo):
        #index in map
        signaldataindex = int(signalinfo["DPOS"])
        samples = self._readSvBy(signaldataindex)
        timesamples = self._getTime(signalinfo)
        return (signalname, samples, timesamples)
           
    #reads signaldata if the layout is FRAME_BY_*, the args are lists
    def _readFrameBy(self, dataposition, repetition, samplecount):
        data = b""
        i = 0
        for frameposition in self.sectionmap["MAP"].values():
            self.fp.seek(frameposition)
            #chunksize = self._readInt()
            self.fp.seek(frameposition + 8 + dataposition)
            if repetition > 1:
                for j in range(samplecount):
                    rawdata = self.fp.read(8)
                    for i in range(repetition):
                        data += rawdata
            else:
                data += self.fp.read(samplecount * 8)
        return np.frombuffer(data, np.float64)

    #reads signaldata if the layout is SV_BY_*
    def _readSvBy(self, signaldataindex):
        data = b''
        self.fp.seek(self.sectionmap["MAP"][signaldataindex])
        datalen = self._readInt()
        if datalen is not None:
            data = self.fp.read(datalen * 8) #assume the dataype is float64
        return np.frombuffer(data, np.float64)

    #generates a time series using sampling time (TC)
    def _genTime(self, tc, datalen):
        data = []
        for i in range(datalen):
            data.append(tc * i)
        return np.array(data)

    #retrieves the signal's time series
    def _getTime(self, signalinfo):
        result = None
        #if the signalinfo contains the info then it acts as a DPOS
        if "TMSV" in signalinfo.keys():
            timesignalindex = int(signalinfo["TMSV"])
            result = self._readSvBy(timesignalindex)
        #if it is in the CONTROL then acts as a signalname
        elif "TMSV" in self.sectionmap["CONTROL"].keys():
            timesignalname = self.sectionmap["CONTROL"]["TMSV"]
            if timesignalname in self.sectionmap["SIGNALS"].keys():
                timesignalindex = self.sectionmap["SIGNALS"][timesignalname]["__INDEX"]
                result = self._readSvBy(timesignalindex)
        #if TC is set in CONTROL
        elif "TC" in self.sectionmap["CONTROL"].keys():
            result = self._genTime(float(self.sectionmap["CONTROL"]["TC"]), int(self.sectionmap["CONTROL"]["NDIM"]))
        if result is None:
            return np.array([])
        else:
            return result

    #reads only one 8 bytes int from the file's actual position
    def _readInt(self):
        result = None
        data = self.fp.read(8)
        if data:
            result = int(struct.unpack('d', data)[0])
        return result
        
    """
    ********* Low level methods end here *******
    ********* Mapping methods start here *******
    """
    #builds a position map about frames or signals
    def _buildMap(self):
        index = 0
        self.fp.seek(self.dataposition)
        numofchunks = self._readInt()
        position = self.dataposition + 8 #data position + number of chunks
        EOF = numofchunks is None
        while not EOF and (index < numofchunks):
            self.fp.seek(position)
            chunksize = self._readInt()
            if chunksize is not None:
                self.sectionmap["MAP"][index] = position
                position += (chunksize + 1) * 8 #assume the datatype is float64, step to the next chunk
                index += 1
            else:
                EOF = True

    #reads the sections
    def _buildSectionmap(self):
        section = ""
        variables = {}
        signalindex = 0
        datareached = False
        EOF = False
        lastsignalname = None
        
        self.fp.seek(0)
        while not datareached and not EOF:
            rawline = self.fp.readline()
            line = rawline.decode('latin-1').rstrip()
            if not rawline:
                #End Of File
                EOF = True
            elif ReadD97.SECTIONEXP.match(line):
                if variables.keys():
                    if section == "CONTROL":
                        self.sectionmap['CONTROL'] = variables
                    elif section == "SIGNAL":
                        if "NAME" in variables.keys():
                            signalname = variables["NAME"]
                            self.sectionmap['SIGNALS'][signalname] = variables
                            self.sectionmap['SIGNALS'][signalname]["__INDEX"] = signalindex
                            self.sectionmap['SIGNALINDEX'][signalindex] = signalname
                            if lastsignalname is not None and "DPOS" in variables.keys():
                                self.sectionmap['SIGNALS'][lastsignalname]["NEXTDPOS"] = variables["DPOS"]
                            lastsignalname = signalname
                    variables = {}
                control = ReadD97.CONTROLEXP.match(line)
                # Is it a signal section?
                signal = ReadD97.SIGNALEXP.match(line)
                # Is it the data section?
                data = ReadD97.DATAEXP.match(line)
                if control is not None:
                    section = "CONTROL"
                elif signal is not None:
                    section = "SIGNAL"
                    signalindex = int(signal.group(1))
                elif data is not None:
                    datareached = True
                else:
                    section = ""
            else:
                #Inside a header section, collect variables
                if section in ["CONTROL", "SIGNAL"]:
                    if line:
                        keyvaluepair = line.split('=')
                        if len(keyvaluepair) > 1:
                            variables[keyvaluepair[0]] = keyvaluepair[1]
                        elif len(keyvaluepair) > 0:
                            variables[keyvaluepair[0]] = None
            self.dataposition = self.fp.tell()
