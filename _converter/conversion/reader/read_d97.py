#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 11:09:05 2019
@author: jez2bp
"""


import re
import struct
import pandas
import numpy
import asammdf


class ReadD97:
    """The class implement functions and methods what can read .D97 files and
    provides the important signal informations for the extracting
    """

    # For section position extracting
    SECTIONEXP = re.compile(r'(\[.*\])')
    CONTROLEXP = re.compile(r'(\[CONTROL\])', re.IGNORECASE)
    SIGNALEXP = re.compile(r'(\[(SIGNAL)(\s*)?\d*\])', re.IGNORECASE)
    DATAEXP = re.compile(r'(\[DATA\])', re.IGNORECASE)

    def __init__(self, filepath):
        self.filepath = filepath
        self.sectionmap = {'CONTROL': {'HEADINDEX': [],
                                       'HEADNAME': [],
                                       'CONTENT': []},
                           'SIGNAL': {'HEADINDEX': [],
                                      'HEADNAME': [],
                                      'CONTENT': [],
                                      'DATA': []},
                           'DATA': {'HEADINDEX': [],
                                    'HEADNAME': [],
                                    'POSITION': 0}}

    def D97Read(self):
        """This is the main method of this class, explore the structure of
        the file, decoded the necessary data and from them provides a proper
        data structure for the conversion.
        """

        # Generating and handling sectionmap variable for make possible the
        # structured data processing
        print('build_sectionmap')
        self.build_sectionmap()
        print('fillwith_defvals')
        self.fillwith_defvals()
        print('build_channeliterators')
        self.build_channeliterators()
        return self.sectionmap

    def get_signal(self, name=None):
        """ This function provides the needed asammdf.MDF() object to the
        Conv4KPI for further processes
        """

        narrowed_mdf = asammdf.MDF()
        if name is None:
            dataframe = self.create_collector_map()
            positions = self.determine_collectible_positions(dataframe)
            collection_dict = self.create_collection_dict(positions)
            signals = self.construct_signal(dataframe, collection_dict)
            print('debug2')
            print(signals)
            for signal in signals:
                narrowed_mdf.append(signal)
        elif isinstance(name, list):
            dataframe = self.create_specific_collector_map(name)
            positions = self.determine_collectible_positions(dataframe)
            collection_dict = self.create_collection_dict(positions)
            signals = self.construct_signal(dataframe, collection_dict)
            for signal in signals:
                narrowed_mdf.append(signal)
        return narrowed_mdf

    def build_sectionmap(self):
        """ The build_sectionmap method reads the file for explore the file
        structure and store necessary information for the successfull parsing.
        """

        # Open the .D97 file to read with buffering
        with open(self.filepath, 'rb', buffering=16384) as file:
            # Initialize some necessary variable
            data_section = False
            control_section = False
            signal_section = False
            temp_dict = {}
            # Read the lines of the file until its reach the data section
            # header line
            for index, rawline in enumerate(file.readlines()):
                if data_section:
                    break
                try:
                    # Decode the lines into latin-1
                    self.sectionmap['DATA']['POSITION'] = self.sectionmap['DATA']['POSITION']+len(rawline)
                    line = rawline.decode('latin-1').rstrip()
                    # Search for string what is matches with any section header
                    # pattern
                    if ReadD97.SECTIONEXP.match(line):
                        # If temp_dict is not empty and any of the switch- like
                        # variables what related to the needed content
                        # structure is active (True) we add the dictionary into
                        # the contenn
                        if temp_dict != {}:
                            if control_section:
                                self.sectionmap['CONTROL']['CONTENT'].append(temp_dict)
                            elif signal_section:
                                self.sectionmap['SIGNAL']['CONTENT'].append(temp_dict)
                        temp_dict = {}
                        # Is it a control section?
                        control = ReadD97.CONTROLEXP.match(line)
                        # Is it a signal section?
                        signal = ReadD97.SIGNALEXP.match(line)
                        # Is it the data section?
                        data = ReadD97.DATAEXP.match(line)
                        if control is not None:
                            self.sectionmap['CONTROL']['HEADINDEX'].append(index)
                            self.sectionmap['CONTROL']['HEADNAME'].append(control.group(1))
                            control_section = True
                            signal_section = False
                        elif signal is not None:
                            self.sectionmap['SIGNAL']['HEADINDEX'].append(index)
                            self.sectionmap['SIGNAL']['HEADNAME'].append(signal.group(1))
                            control_section = False
                            signal_section = True
                        elif data is not None:
                            self.sectionmap['DATA']['HEADINDEX'].append(index)
                            self.sectionmap['DATA']['HEADNAME'].append(data.group(1))
                            data_section = True
                            control_section = False
                            signal_section = False
                        else:
                            control_section = False
                            signal_section = False
                    # If the examined line is not a section header but the
                    # previos section header was one of the interesting
                    # sections we split the key-value pairs into parts and
                    # store it in the temp_dict
                    else:
                        if signal_section or control_section:
                            kvp = line.split('=')
                            if kvp[0]:
                                try:
                                    temp_dict[kvp[0]] = kvp[1]
                                except IndexError:
                                    temp_dict[kvp[0]] = None
                except UnicodeDecodeError:
                    pass

    def build_channeliterators(self):
        """
        The function generates the iterator object of the signal samples
        """

        with open(self.filepath, 'rb', buffering=16384) as file:
            file.seek(int(self.sectionmap['DATA']['POSITION']))
            segment = file.read(8)
            n = int(struct.unpack('d', segment)[0])
            file.seek(int(self.sectionmap['DATA']['POSITION'])+8)
            m_n_list = []
            m_n_pos = []
            channel_iters = []
            ok = True
            dpos = 0
            bpos = 0
            while ok:
                r = file.read(8)
                bpos = bpos + len(r)
                ok = r != b'' and dpos < n
                if ok:
                    m_n = int(struct.unpack('d', r)[0])
                    m_n_list.append(m_n)
                    file.seek(m_n*8, 1)
                    m_n_pos.append(bpos+int(self.sectionmap['DATA']['POSITION'])+8)
                    bpos = bpos + m_n*8
                dpos = dpos+1
            for index, byte_pos in enumerate(m_n_pos):
                file.seek(byte_pos)
                raw = file.read(8*m_n_list[index])
                channel_iter = struct.iter_unpack('d', raw)
                channel_iters.append(channel_iter)
        self.sectionmap['SIGNAL']['DATA'] = channel_iters

    def fillwith_defvals(self):
        """The functions loads some default necessary information into the
        sectionmap if its not available
        """
        # CONTROL SECTION
        c_TC = 0
        c_TMSV = None
        if 'TC' not in self.sectionmap['CONTROL']['CONTENT'][0].keys():
            self.sectionmap['CONTROL']['CONTENT'][0]['TC'] = c_TC
        if 'TMSV' not in self.sectionmap['CONTROL']['CONTENT'][0].keys():
            self.sectionmap['CONTROL']['CONTENT'][0]['TMSV'] = c_TMSV
        # SIGNAL SECTION
        s_RPT = 1
        for index, signal_cont in enumerate(self.sectionmap['SIGNAL']['CONTENT']):
            if 'RPT' not in signal_cont.keys():
                self.sectionmap['SIGNAL']['CONTENT'][index]['RPT'] = s_RPT

    def create_collector_map(self):
        """The function creates a dataframe with the necessary information of
        sample data extracting
        """

        indexes = []
        names = []
        TIMEPOS = []
        DPOS = []
        for index, item in enumerate(self.sectionmap['SIGNAL']['CONTENT']):
            if 'TMSV' in item.keys():
                indexes.append(index)
                names.append(item['NAME'])
                TIMEPOS.append(int(item['TMSV']))
                if 'DPOS' in item.keys():
                    DPOS.append(int(item['DPOS']))
                else:
                    DPOS.append(None)
            elif 'TMSV' in self.sectionmap['CONTROL']['CONTENT'][0].keys():
                if self.sectionmap['CONTROL']['CONTENT'][0]['TMSV'] is not None:
                    indexes.append(int(index))
                    names.append(item['NAME'])
                    for index, signaldatas in enumerate(self.sectionmap['SIGNAL']['CONTENT']):
                        signalname = self.sectionmap['SIGNAL']['CONTENT'][index]['NAME']
                        if signalname == self.sectionmap['CONTROL']['CONTENT'][0]['TMSV']:
                            TIMEPOS.append(int(self.sectionmap['SIGNAL']['CONTENT'][index]['DPOS']))
                            break
                    if 'DPOS' in item.keys():
                        DPOS.append(int(item['DPOS']))
                    else:
                        DPOS.append(None)
                else:#TODO quick dirty solution for missing tmsv
                    indexes.append(index)
                    names.append(item['NAME'])
                    TIMEPOS.append(int(0))
                    if 'DPOS' in item.keys():
                        DPOS.append(int(item['DPOS']))
                    else:
                        DPOS.append(None)
                
        dataframe = pandas.DataFrame({'INDEX':indexes,
                                      'NAME':names,
                                      'TIMELOC': TIMEPOS,
                                      'DATALOC': DPOS})
        return dataframe

    def create_specific_collector_map(self, namelist):
        """The function creates a dataframe with the necessary information of
        sample data extracting but just with the needed signals
        """

        namelist = list(set(namelist))
        indexes = []
        names = []
        TIMEPOS = []
        DPOS = []
        for name in namelist:
            for index, item in enumerate(self.sectionmap['SIGNAL']['CONTENT']):
                if item['NAME'] == name:
                    if 'TMSV' in item.keys():
                        indexes.append(index)
                        names.append(item['NAME'])
                        TIMEPOS.append(int(item['TMSV']))
                        if 'DPOS' in item.keys():
                            DPOS.append(int(item['DPOS']))
                        else:
                            DPOS.append(None)
                    elif 'TMSV' in self.sectionmap['CONTROL']['CONTENT'][0].keys():
                        if self.sectionmap['CONTROL']['CONTENT'][0]['TMSV'] is not None:
                            indexes.append(index)
                            names.append(item['NAME'])
                            for index, signaldatas in enumerate(self.sectionmap['SIGNAL']['CONTENT']):
                                signalname = self.sectionmap['SIGNAL']['CONTENT'][index]['NAME']
                                if signalname == self.sectionmap['CONTROL']['CONTENT'][0]['TMSV']:
                                    TIMEPOS.append(int(self.sectionmap['SIGNAL']['CONTENT'][index]['DPOS']))
                                    break
                            if 'DPOS' in item.keys():
                                DPOS.append(int(item['DPOS']))
                            else:
                                DPOS.append(None)
                        else:#TODO quick dirty solution for missing tmsv
                            indexes.append(index)
                            names.append(item['NAME'])
                            TIMEPOS.append(int(0))
                            if 'DPOS' in item.keys():
                                DPOS.append(int(item['DPOS']))
                            else:
                                DPOS.append(None)
        dataframe = pandas.DataFrame({'INDEX':indexes,
                                      'NAME':names,
                                      'TIMELOC': TIMEPOS,
                                      'DATALOC': DPOS})
        return dataframe

    def determine_collectible_positions(self, dataframe):
        '''The function determines all postions what need to be collected'''

        positions_to_collect = []
        for pos in dataframe['TIMELOC']:
            positions_to_collect.append(int(pos))
        for pos in dataframe['DATALOC']:
            positions_to_collect.append(int(pos))
        positions_to_collect = list(set(positions_to_collect))
        return positions_to_collect

    def create_collection_dict(self, positions_to_collect):
        """The function creates a dictionary with the necessary data for the
        signal construction
        """

        collection_dict = {}
        for pos in positions_to_collect:
            collection_dict[pos] = []
            iterator = self.sectionmap['SIGNAL']['DATA'][pos]
            for value in iterator:
                collection_dict[pos].append(value[0])
        return collection_dict

    def construct_signal(self, dataframe, collection_dict):
        """The function creates a list of asammdf signal, each of the signals
        will be in a different signal group
        """

        signals = []
        for index, signal_name in enumerate(dataframe['NAME']):
            time = numpy.array(collection_dict[int(dataframe['TIMELOC'][index])],
                               dtype=float)
            sample = numpy.array(collection_dict[int(dataframe['DATALOC'][index])],
                                 dtype=float)
            name = signal_name
            signal = asammdf.mdf.Signal(name=name,
                                        samples=sample,
                                        timestamps=time)
            signals.append([signal])
        return signals
