#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 12:42:08 2019

@author: jez2bp
"""

import conversion.reader.read_d97
import struct
import unittest
import asammdf
import numpy
import pandas

class TestReadD97(unittest.TestCase):

    def setUp(self):

        self.filepath = "tests/testD97.d97"
        self.signalnames = ['A', 'B', 'C']

        self.code = conversion.reader.read_d97.ReadD97(self.filepath)

    def test_D97Read(self):
        my_sectionmap = self.code.D97Read()
        self.assertIsInstance(my_sectionmap, dict)
        self.assertTrue('CONTROL' in my_sectionmap.keys() and 'DATA' in my_sectionmap.keys() and 'SIGNAL' in my_sectionmap.keys())

        self.assertIsInstance(my_sectionmap['CONTROL'], dict)
        self.assertIsInstance(my_sectionmap['DATA'], dict)
        self.assertIsInstance(my_sectionmap['SIGNAL'], dict)
        self.assertTrue('HEADINDEX' in my_sectionmap['CONTROL'].keys() and 'HEADNAME' in my_sectionmap['CONTROL'].keys() and 'CONTENT' in my_sectionmap['CONTROL'].keys())
        self.assertTrue('HEADINDEX' in my_sectionmap['SIGNAL'].keys() and 'HEADNAME' in my_sectionmap['SIGNAL'].keys() and 'CONTENT' in my_sectionmap['SIGNAL'].keys() and 'DATA' in my_sectionmap['SIGNAL'].keys())
        self.assertTrue('HEADINDEX' in my_sectionmap['DATA'].keys() and 'HEADNAME' in my_sectionmap['DATA'].keys() and 'POSITION' in my_sectionmap['DATA'].keys())

        self.assertIsInstance(my_sectionmap['CONTROL']['HEADINDEX'], list)
        self.assertIsInstance(my_sectionmap['CONTROL']['HEADNAME'], list)
        self.assertIsInstance(my_sectionmap['CONTROL']['CONTENT'], list)
        self.assertIsInstance(my_sectionmap['SIGNAL']['HEADINDEX'], list)
        self.assertIsInstance(my_sectionmap['SIGNAL']['HEADNAME'], list)
        self.assertIsInstance(my_sectionmap['SIGNAL']['CONTENT'], list)
        self.assertIsInstance(my_sectionmap['SIGNAL']['DATA'], list)
        self.assertIsInstance(my_sectionmap['DATA']['HEADINDEX'], list)
        self.assertIsInstance(my_sectionmap['DATA']['HEADNAME'], list)
        self.assertIsInstance(my_sectionmap['DATA']['POSITION'], int)
        self.assertTrue(len(my_sectionmap['CONTROL']['HEADINDEX'])==1)
        self.assertIsInstance(my_sectionmap['CONTROL']['HEADINDEX'][0], int)
        self.assertTrue(len(my_sectionmap['CONTROL']['HEADINDEX']) == len(my_sectionmap['CONTROL']['HEADNAME']))
        self.assertEqual('[CONTROL]', my_sectionmap['CONTROL']['HEADNAME'][0])
        self.assertTrue(len(my_sectionmap['CONTROL']['HEADNAME']) == len(my_sectionmap['CONTROL']['CONTENT']))
        self.assertIsInstance(my_sectionmap['CONTROL']['CONTENT'][0], dict)

        self.assertTrue(len(my_sectionmap['SIGNAL']['HEADINDEX']) == len(my_sectionmap['SIGNAL']['HEADNAME']) == len(my_sectionmap['SIGNAL']['CONTENT']) == len(my_sectionmap['SIGNAL']['DATA']))

        hitype = []
        for element in my_sectionmap['SIGNAL']['HEADINDEX']:
            hitype.append(type(element))
        self.assertTrue(len(list(set(hitype))), 1)
        self.assertIsInstance(my_sectionmap['SIGNAL']['HEADINDEX'][0], int)

        hntype = []
        for element in my_sectionmap['SIGNAL']['HEADNAME']:
            hntype.append(type(element))
        self.assertTrue(len(list(set(hntype))), 1)
        self.assertIsInstance(my_sectionmap['SIGNAL']['HEADNAME'][0], str)

        conttype = []
        for element in my_sectionmap['SIGNAL']['CONTENT']:
            conttype.append(type(element))
        self.assertTrue(len(list(set(conttype))), 1)
        self.assertIsInstance(my_sectionmap['SIGNAL']['CONTENT'][0], dict)

        datatype = []
        for element in my_sectionmap['SIGNAL']['DATA']:
            datatype.append(type(element))
        self.assertTrue(len(list(set(datatype))), 1)
        self.assertIsInstance(my_sectionmap['SIGNAL']['DATA'][0], type(struct.iter_unpack('d', b'abcabcab')))

        self.assertTrue(len(my_sectionmap['DATA']['HEADINDEX']) == len(my_sectionmap['DATA']['HEADNAME']) == 1)
        self.assertIsInstance(my_sectionmap['DATA']['HEADINDEX'][0], int)
        self.assertEqual(my_sectionmap['DATA']['HEADNAME'][0], '[DATA]')
        self.assertIsInstance(my_sectionmap['DATA']['POSITION'], int)

    def test_get_signal(self):
        self.code.D97Read()

        my_mdf = self.code.get_signal(self.signalnames)
        self.assertIsInstance(my_mdf,
                              type(asammdf.mdf.MDF()))

        chaniter = my_mdf.iter_channels()
        names = []
        samples = []
        times = []
        for index, sig in enumerate(chaniter):
            names.append(sig.name)
            samples.append(sig.samples)
            times.append(sig.timestamps)

        self.assertTrue(len(names) == len(self.signalnames) == len(samples) == len(times) == 3)
        self.assertEqual(sorted(self.signalnames),
                         sorted(names))

        samptype = []
        for element in samples:
            samptype.append(type(element))
        self.assertTrue(len(list(set(samptype)))==1)
        self.assertIsInstance(samples[0], type(numpy.array([0])))

        timetype = []
        for element in times:
            timetype.append(type(element))
        self.assertTrue(len(list(set(timetype)))==1)
        self.assertIsInstance(times[0], type(numpy.array([0])))

    def test_create_collector_map(self):

        self.code.D97Read()
        my_frame = self.code.create_collector_map()
        self.assertIsInstance(my_frame, type(pandas.DataFrame()))
        self.assertTrue('INDEX' in list(my_frame) and 'NAME' in list(my_frame) and 'TIMELOC' in list(my_frame) and 'DATALOC' in list(my_frame))

        names = []
        timelocs = []
        datalocs = []
        for index, item in enumerate(my_frame['INDEX']):
            names.append(my_frame['NAME'][index])
            timelocs.append(my_frame['TIMELOC'][index])
            datalocs.append(my_frame['DATALOC'][index])
        self.assertTrue('A' in names and 'B' in names and 'C' in names and 'D' in names and 'E' in names)
        self.assertTrue('TIME' in names) #the test .d97 generated by mdfdset6c from the testMDF.mf4 and the conversion extends the signal list with the TIME signal which is the common TIME data of all signals

        timetypes = []
        for element in timelocs:
            timetypes.append(type(element))
        self.assertTrue(len(list(set(timetypes))) == 1)
        self.assertIsInstance(timelocs[0], numpy.int64)

        datatypes = []
        for element in datalocs:
            datatypes.append(type(element))
        self.assertTrue(len(list(set(datatypes))) == 1)
        self.assertIsInstance(datalocs[0], numpy.int64)

    def test_create_specific_collector_map(self):

        self.code.D97Read()
        my_frame = self.code.create_specific_collector_map(self.signalnames)
        self.assertIsInstance(my_frame, type(pandas.DataFrame()))
        self.assertTrue('INDEX' in list(my_frame) and 'NAME' in list(my_frame) and 'TIMELOC' in list(my_frame) and 'DATALOC' in list(my_frame))
        names = []
        timelocs = []
        datalocs = []
        for index, item in enumerate(my_frame['INDEX']):
            names.append(my_frame['NAME'][index])
            timelocs.append(my_frame['TIMELOC'][index])
            datalocs.append(my_frame['DATALOC'][index])
        self.assertTrue('A' in names and 'B' in names and 'C' in names and 'D' not in names and 'E' not in names)
        self.assertTrue('TIME' not in names) #the test .d97 generated by mdfdset6c from the testMDF.mf4 and the conversion extends the signal list with the TIME signal which is the common TIME data of all signals

        timetypes = []
        for element in timelocs:
            timetypes.append(type(element))
        self.assertTrue(len(list(set(timetypes))) == 1)
        self.assertIsInstance(timelocs[0], numpy.int64)

        datatypes = []
        for element in datalocs:
            datatypes.append(type(element))
        self.assertTrue(len(list(set(datatypes))) == 1)
        self.assertIsInstance(datalocs[0], numpy.int64)

    def test_determine_collectible_positions(self):
        data = {'INDEX': [1,2,3],
                'NAME': ['A', 'B', 'C'],
                'TIMELOC': [1,2,3],
                'DATALOC': [4,5,6]}
        my_frame = pandas.DataFrame(data)
        my_positions = self.code.determine_collectible_positions(my_frame)
        self.assertIsInstance(my_positions, list)
        postypes = []
        for element in my_positions:
            postypes.append(type(element))
        self.assertTrue(len(list(set(postypes))) == 1)
        self.assertIsInstance(my_positions[0], int)
        self.assertEqual([1,2,3,4,5,6], sorted(my_positions))

    def test_create_collection_dict(self):
        self.code.D97Read()
        my_positions = [0,1,2,3,4,5]
        my_collection = self.code.create_collection_dict(my_positions)
        self.assertIsInstance(my_collection, dict)
        self.assertTrue(0 in my_collection.keys() and 1 in my_collection.keys() and 2 in my_collection.keys() and 3 in my_collection.keys() and 4 in my_collection.keys() and 5 in my_collection.keys())
        self.assertIsInstance(my_collection[0], list)
        self.assertIsInstance(my_collection[1], list)
        self.assertIsInstance(my_collection[2], list)
        self.assertIsInstance(my_collection[3], list)
        self.assertIsInstance(my_collection[4], list)
        self.assertIsInstance(my_collection[5], list)

        mc0types = []
        for element in my_collection[0]:
            mc0types.append(type(element))
        self.assertTrue(len(list(set(mc0types))) == 1)
        self.assertIsInstance(my_collection[0][0], float)

        mc1types = []
        for element in my_collection[1]:
            mc1types.append(type(element))
        self.assertTrue(len(list(set(mc1types))) == 1)
        self.assertIsInstance(my_collection[1][0], float)

        mc2types = []
        for element in my_collection[2]:
            mc2types.append(type(element))
        self.assertTrue(len(list(set(mc2types))) == 1)
        self.assertIsInstance(my_collection[2][0], float)

        mc3types = []
        for element in my_collection[3]:
            mc3types.append(type(element))
        self.assertTrue(len(list(set(mc3types))) == 1)
        self.assertIsInstance(my_collection[3][0], float)

        mc4types = []
        for element in my_collection[4]:
            mc4types.append(type(element))
        self.assertTrue(len(list(set(mc4types))) == 1)
        self.assertIsInstance(my_collection[4][0], float)

        mc5types = []
        for element in my_collection[5]:
            mc5types.append(type(element))
        self.assertTrue(len(list(set(mc5types))) == 1)
        self.assertIsInstance(my_collection[5][0], float)

    def test_construct_signal(self):
        data = {'INDEX':[1,2,3],
                'NAME':['A', 'B', 'C'],
                'TIMELOC': [0,1,2],
                'DATALOC':[3,4,5]}
        my_frame = pandas.DataFrame(data)
        my_collection = {0:numpy.arange(0,10,1),
                         1:numpy.arange(0,10,2),
                         2:numpy.arange(0,30,1),
                         3:numpy.arange(0,10,1),
                         4:numpy.arange(0,10,2),
                         5:numpy.arange(0,30,1)}
        my_signals = self.code.construct_signal(my_frame, my_collection)
        self.assertIsInstance(my_signals, list)
        mstypes = []
        for element in my_signals:
            mstypes.append(type(element))
        self.assertTrue(len(list(set(mstypes)))==1)
        self.assertIsInstance(my_signals[0], list)
        self.assertIsInstance(my_signals[1], list)
        self.assertIsInstance(my_signals[2], list)

        self.assertTrue(len(my_signals[0])==len(my_signals[1])==len(my_signals[2])==1)

        refsamples = numpy.arange(0,2,1)
        reftimes = numpy.arange(0,2,1)
        refname = 'REFNAME'

        self.assertIsInstance(my_signals[0][0], type(asammdf.mdf.Signal(name = refname, timestamps = reftimes, samples = refsamples)))
        self.assertIsInstance(my_signals[1][0], type(asammdf.mdf.Signal(name = refname, timestamps = reftimes, samples = refsamples)))
        self.assertIsInstance(my_signals[2][0], type(asammdf.mdf.Signal(name = refname, timestamps = reftimes, samples = refsamples)))

        self.assertEqual(my_signals[0][0].name, 'A')
        self.assertEqual(my_signals[1][0].name, 'B')
        self.assertEqual(my_signals[2][0].name, 'C')

        self.assertTrue(list(my_signals[0][0].samples) == list(my_signals[0][0].timestamps) == list(numpy.arange(0,10,1)))
        self.assertTrue(list(my_signals[1][0].samples) == list(my_signals[1][0].timestamps) == list(numpy.arange(0,10,2)))
        self.assertTrue(list(my_signals[2][0].samples) == list(my_signals[2][0].timestamps) == list(numpy.arange(0,30,1)))
