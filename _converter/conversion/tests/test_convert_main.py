#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 15:34:53 2019

@author: jez2bp
"""
import os
import scipy.io
import numpy
import unittest
import conversion.convert_main
import scipy
import asammdf
import threading

class TestConvertMain(unittest.TestCase):

    def __init__(self, file, force, evalsignals, testName):
        super(TestConvertMain, self).__init__(testName)
        self.file = file
        self.force = force
        self.evalsignals = evalsignals

    def setUp(self):
        self.pltpath = 'tests/testPLT.PLT'
        self.inputfilename = self.file
        self.outputfilename = 'tests/convertmainTESTMAT.mat'
        self.m_eval = "testeval"
        threadswitch = threading.Event()
        threadswitch.clear()
        gui = False
        self.code = conversion.convert_main.MainConversion(self.pltpath,
                                                           self.inputfilename,
                                                           self.outputfilename,
                                                           self.force,
                                                           self.evalsignals,
                                                           self.m_eval,
                                                           threadswitch,
                                                           gui,
                                                           None,
                                                           None)

    def test_start_conversion(self):
        res = self.code.start_conversion()

        self.assertIsInstance(res, list, "MainConversion.start_conversion: checking rerturn value type")

        name = res[0]

        self.assertIsInstance(name, str, "MainConversion.start_conversion: checking rerturn list element type")
        self.assertEqual(name, 'tests/convertmainTESTMAT.mat', "MainConversion.start_conversion: checking rerturn list element value")
        self.assertTrue(os.path.isfile(name), "MainConversion.start_conversion: checking if the file was created successfully")

        mymat = scipy.io.loadmat(name)
        myref = self.dictionary_data_for_testing()
        signalnameresults = []
        signalsampleresults = []
        lensofsamples = []
        for key in myref.keys():
            signalnameresults.append(key in mymat.keys())
            signalsampleresults.append(mymat[key][0])
            lensofsamples.append(len(mymat[key][0]))

        self.assertTrue(list(set(signalnameresults))[0], "MainConversion.start_conversion: checking if the file contains all signalname")
        self.assertTrue(len(list(set(lensofsamples))) == 1 and list(set(lensofsamples))[0] != 0, "MainConversion.start_conversion: checking if all sample has the same length")

        compare = []
        for sample in signalsampleresults:
            for val in sample:
                compare.append(isinstance(val, int) or isinstance(val, float))

        self.assertTrue(list(set(compare))[0], "MainConversion.start_conversion: checking if all sample value has the right type")
        os.remove(name)

    def test_change_names_and_interpolate(self):
        sourcenames = []
        ref_mdf = self.reference_mdf()
        refitered = ref_mdf.iter_channels(skip_master=True)
        for signal in refitered:
            sourcenames.append(signal.name)
        pltdict = self.dictionary_data_for_testing()

        self.assertIsInstance(self.code.change_names_and_interpolate(ref_mdf, pltdict),
                              type(asammdf.mdf.MDF()),
                              "MainConversion.change_names_and_interpolate: checking return value type")

        my_mdf = self.code.change_names_and_interpolate(ref_mdf, pltdict)
        targetnames = []
        samplens = []
        stamplens = []
        targitered = my_mdf.iter_channels(skip_master=True)
        for signal in targitered:
            targetnames.append(signal.name)
            samplens.append(len(list(signal.samples)))
            stamplens.append(len(list(signal.timestamps)))

        self.assertEqual(sorted(sourcenames),
                         sorted(['B', 'C', 'D', 'E', 'A']),
                         "MainConversion.change_names_and_interpolate: checking if the reference mdf contains all sourcenames")
        self.assertEqual(sorted(['A_1', 'A_2', 'B_3', 'C_4', 'D_5', 'E_6', 'B_7', 'D_8', 'E_9']),
                         sorted(targetnames),
                         "MainConversion.change_names_and_interpolate: checking if the renamed mdf contains all tergetnames after the process")
        self.assertEqual(len(list(set(samplens))),
                         len(list(set(stamplens))),
                         "MainConversion.change_names_and_interpolate: checking if all sample has the same length after the process")

    def test_generate_new_signals_and_rename(self):
        sourcenames = []
        ref_mdf = self.reference_mdf()
        refitered = ref_mdf.iter_channels(skip_master=True)
        for signal in refitered:
            sourcenames.append(signal.name)
        pltdict = self.dictionary_data_for_testing()

        self.assertIsInstance(self.code.generate_new_signals_interp_and_rename(ref_mdf,
                                                                               pltdict),
                              type(asammdf.mdf.MDF()),
                              "MainConversion.generate_new_signals_and_rename: checking return value type")

        my_mdf = self.code.generate_new_signals_interp_and_rename(ref_mdf, pltdict)
        targetnames = []
        samplens = []
        stamplens = []
        targitered = my_mdf.iter_channels(skip_master=True)
        for signal in targitered:
            targetnames.append(signal.name)
            samplens.append(len(list(signal.samples)))
            stamplens.append(len(list(signal.timestamps)))

        self.assertEqual(sorted(sourcenames),
                         sorted(['B', 'C', 'D', 'E', 'A']),
                         "MainConversion.generate_new_signals_and_rename: checking if the renamed mdf contains all tergetnames after the process")
        self.assertEqual(sorted(['A_1', 'A_2', 'B_3', 'C_4', 'D_5', 'E_6', 'B_7', 'D_8', 'E_9']),
                         sorted(targetnames),
                         "MainConversion.generate_new_signals_and_rename: checking if the renamed mdf contains all tergetnames after the process")
        self.assertEqual(len(list(set(samplens))),
                         len(list(set(stamplens))),
                         "MainConversion.generate_new_signals_and_rename: checking if all sample has the same length after the process")

    def test_get_plt_data(self):
        my_pltdict, my_sourcenames = self.code.get_plt_data()

        self.assertIsInstance(my_pltdict, dict)
        self.assertIsInstance(my_sourcenames, list)

        reference_dict = self.dictionary_data_for_testing()

        self.assertEqual(my_pltdict, reference_dict)
        self.assertEqual(sorted(list(set(my_sourcenames))), ['A', 'B', 'C', 'D', 'E'])

    def test_get_data(self):
        my_sourcenames = ['A', 'A', 'B', 'C', 'D', 'E', 'B', 'D', 'E']
        my_data = self.code.get_data(['A', 'B', 'C', 'D', 'E'],
                                     self.dictionary_data_for_testing(),
                                     self.inputfilename)

        self.assertIsInstance(my_data, type(asammdf.MDF()))

        datatypes = []
        for element in my_data:
            datatypes.append(type(element))

        self.assertTrue(len(list(set(datatypes))) == 1)
        self.assertIsInstance(my_data, type(asammdf.MDF()))

        mdf_signals = []
        signaltypes = []
        chaniter = my_data.iter_channels(skip_master=True)
        for signal in chaniter:
            mdf_signals.append(signal)
            signaltypes.append(type(signal))

        self.assertTrue(len(list(set(signaltypes))) == 1)
        self.assertIsInstance(mdf_signals[0],
                              type(asammdf.Signal(name="REFNAME",
                                                  samples=numpy.arange(0, 10, 1),
                                                  timestamps=numpy.arange(0, 10, 1))))
        self.assertEqual(len(list(set(my_sourcenames))),
                         len(mdf_signals))

    def test_remove_extension(self):
        self.assertEqual(self.code.remove_extension('Z:\CC\CC-Bp\Public'),
                         'Z:\CC\CC-Bp\Public')
        self.assertEqual(self.code.remove_extension('Z:\CC\CC-Bp\Public\mymeas'),
                         'Z:\CC\CC-Bp\Public\mymeas')
        self.assertEqual(self.code.remove_extension('Z:\CC\CC-Bp\Public\parameters.cool.json'),
                         'Z:\CC\CC-Bp\Public\parameters.cool.json')
        self.assertEqual(self.code.remove_extension('Z:\CC\CC-Bp\Public\parking.example.mf4'),
                         'Z:\CC\CC-Bp\Public\parking.example')
        self.assertEqual(self.code.remove_extension('Z:\CC\CC-Bp\Public\parking.example.d97'),
                         'Z:\CC\CC-Bp\Public\parking.example')
        self.assertEqual(self.code.remove_extension('Z:\CC\CC-Bp\Public\parking.example.dat'),
                         'Z:\CC\CC-Bp\Public\parking.example')
        self.assertEqual(self.code.remove_extension('Z:\CC\CC-Bp\Public\parking.example.MdF'),
                         'Z:\CC\CC-Bp\Public\parking.example')
        self.assertEqual(self.code.remove_extension('Z:\CC\CC-Bp\Public\parking.example.mat'),
                         'Z:\CC\CC-Bp\Public\parking.example.mat')

    def test_get_dat_mdf_data(self):
        if self.code.is_extension(self.file, ['dat', 'mf4', 'mdf']):
            my_sourcenames = ['A', 'A', 'B', 'C', 'D', 'E', 'B', 'D', 'E']
            pltdict, sourcenames = self.code.get_plt_data()
            getdata_result = self.code.get_dat_mdf_data(self.file, sourcenames, pltdict)

            self.assertIsInstance(getdata_result, dict)
            self.assertEqual(self.file, getdata_result['filename'])
            self.assertTrue(getdata_result['success'])
            self.assertFalse(getdata_result['stopped'])
            self.assertFalse(getdata_result['skip'])
            self.assertTrue(getdata_result['error_msg'] is None)
            self.assertEqual(sorted(getdata_result['mdf_channels']), sorted(list(set(sourcenames))))

            my_data = getdata_result['mdf_object']

            self.assertIsInstance(my_data, type(asammdf.MDF()))

            datatypes = []
            for element in my_data:
                datatypes.append(type(element))

            self.assertTrue(len(list(set(datatypes))) == 1)
            self.assertIsInstance(my_data, type(asammdf.MDF()))

            mdf_signals = []
            signaltypes = []
            chaniter = my_data.iter_channels(skip_master=True)
            for signal in chaniter:
                mdf_signals.append(signal)
                signaltypes.append(type(signal))

            self.assertTrue(len(list(set(signaltypes))) == 1)
            self.assertIsInstance(mdf_signals[0],
                                  type(asammdf.Signal(name="REFNAME",
                                                      samples=numpy.arange(0, 10, 1),
                                                      timestamps=numpy.arange(0, 10, 1))))
            self.assertEqual(len(list(set(my_sourcenames))),
                             len(mdf_signals))

    def test_get_d97_data(self):
        if self.code.is_extension(self.file, ['d97']):
            my_sourcenames = ['A', 'A', 'B', 'C', 'D', 'E', 'B', 'D', 'E']
            pltdict, sourcenames = self.code.get_plt_data()
            getdata_result = self.code.get_d97_data(self.file,
                                                    sourcenames,
                                                    pltdict)

            self.assertIsInstance(getdata_result, dict)
            self.assertEqual(self.file, getdata_result['filename'])
            self.assertTrue(getdata_result['success'])
            self.assertFalse(getdata_result['stopped'])
            self.assertFalse(getdata_result['skip'])
            self.assertTrue(getdata_result['error_msg'] is None)
            self.assertEqual(sorted(getdata_result['mdf_channels']),
                             sorted(list(set(sourcenames))))

            my_data = getdata_result['mdf_object']

            self.assertIsInstance(my_data, type(asammdf.MDF()))

            datatypes = []
            for element in my_data:
                datatypes.append(type(element))

            self.assertTrue(len(list(set(datatypes))) == 1)
            self.assertIsInstance(my_data, type(asammdf.MDF()))

            mdf_signals = []
            signaltypes = []
            chaniter = my_data.iter_channels(skip_master=True)
            for signal in chaniter:
                mdf_signals.append(signal)
                signaltypes.append(type(signal))

            self.assertTrue(len(list(set(signaltypes))) == 1)
            self.assertIsInstance(mdf_signals[0],
                                  type(asammdf.Signal(name="REFNAME",
                                                      samples=numpy.arange(0, 10, 1),
                                                      timestamps=numpy.arange(0, 10, 1))))
            self.assertEqual(len(list(set(my_sourcenames))),
                             len(mdf_signals))

    def reference_mdf(self):
        names = ['B', 'C', 'D', 'E', 'A']
        timestamps = [numpy.arange(0, 10, 2),
                      numpy.arange(0, 10, 1),
                      numpy.arange(0, 10, 3),
                      numpy.arange(5, 20, 1),
                      numpy.arange(-100, 8, 1)]
        samples = [numpy.arange(0, 10, 2),
                   numpy.arange(0, 10, 1),
                   numpy.arange(0, 10, 3),
                   numpy.arange(5, 20, 1),
                   numpy.arange(-100, 8, 1)]
        ref_mdf = asammdf.mdf.MDF(memory='full')
        for index, name in enumerate(names):
            signal = [asammdf.mdf.Signal(samples=samples[index],
                                         name=name,
                                         timestamps=timestamps[index])]
            ref_mdf.append(signal)
        return ref_mdf

    def dictionary_data_for_testing(self):
        test_dict = {'A_1': {'scmin': -1000.0,
                             'scmax': 1000.0,
                             'name': 'A',
                             'unit': 'potatoe',
                             'bit': 0,
                             'nobit': 1,
                             'exp': ''},
                     'A_2': {'scmin': -1000.0,
                             'scmax': 1000.0,
                             'name': 'A',
                             'unit': 'tomatoe',
                             'bit': 1,
                             'nobit': 1,
                             'exp': '+ exp=x/2'},
                     'B_3': {'scmin': -800.0,
                             'scmax': 800.0,
                             'name': 'B',
                             'unit': 'suffer',
                             'bit': 0,
                             'nobit': 1,
                             'exp': '+ exp(x1=A_2)=(x/2)+x1/x-1'},
                     'C_4': {'scmin': -200.0,
                             'scmax': 200.0,
                             'name': 'C',
                             'unit': 'pain',
                             'bit': 1,
                             'nobit': 1,
                             'exp': ''},
                     'D_5': {'scmin': -150.0,
                             'scmax': 150.0,
                             'name': 'D',
                             'unit': 'pony',
                             'bit': 0,
                             'nobit': 1,
                             'exp': '+ exp(x1=B_3,x2=A_2)=(x1+x2+x)/3'},
                     'E_6': {'scmin': -300.0,
                             'scmax': 300.0,
                             'name': 'E',
                             'unit': 'dragon',
                             'bit': 1,
                             'nobit': 2,
                             'exp': ''},
                     'B_7': {'scmin': -250.0,
                             'scmax': 250.0,
                             'name': 'B',
                             'unit': 'the_winter_is_coming',
                             'bit': 0,
                             'nobit': 5,
                             'exp': ''},
                     'D_8': {'scmin': -1000.0,
                             'scmax': 1000.0,
                             'name': 'D',
                             'unit': 'a_lannister_always_pays_his_debts',
                             'bit': 1,
                             'nobit': 8,
                             'exp': '+ exp=x/3'},
                     'E_9': {'scmin': -70.0,
                             'scmax': 7.0,
                             'name': 'E',
                             'unit': 'for_frodo',
                             'bit': 0,
                             'nobit': 100,
                             'exp': ''}}
        return test_dict
