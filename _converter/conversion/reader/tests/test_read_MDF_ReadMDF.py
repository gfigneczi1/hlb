#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 16:29:03 2019

@author: jez2bp
"""
import numpy
import asammdf
import unittest
import conversion.reader.read_mdf


class TestReadMDF(unittest.TestCase):

    def __init__(self, testName):
        super(TestReadMDF, self).__init__(testName)

    def setUp(self):

        self.filepath = "tests/testMDF.mf4"
        self.plt_signals = ['A', 'B', 'C', 'D', 'E']

        self.code = conversion.reader.read_mdf.ReadMDF(filepath=self.filepath,
                                                       plt_signals=self.plt_signals)

    def test_read_channels(self):

        selected_channels = ['A', 'B', 'C', 'D', 'E', 'A']
        self.assertIsInstance(self.code.read_channels(),
                              tuple)
        self.assertIsInstance(self.code.read_channels()[0],
                              type(asammdf.mdf.MDF()))
        self.assertIsInstance(self.code.read_channels()[1],
                              list)

        my_mdf = self.code.read_channels()[0]
        names_from_mdf, times_from_mdf, samples_from_mdf = self.generate_comperable_data_from_read(my_mdf)
        test_channels, timestamps_for_test, samples_for_test = self.reference_data()

        self.assertEqual(sorted(names_from_mdf), sorted(test_channels))
        self.assertEqual(sorted(times_from_mdf), sorted(timestamps_for_test))
        self.assertEqual(sorted(samples_from_mdf), sorted(samples_for_test))
        self.assertEqual(sorted(list(set(selected_channels))), sorted(test_channels))
        self.assertEqual(sorted(self.code.read_channels()[1]), sorted(test_channels))
        selected_channels.append('X')
        my_mdf = self.code.read_channels()[0]
        names_from_mdf, times_from_mdf, samples_from_mdf = self.generate_comperable_data_from_read(my_mdf)
        test_channels, timestamps_for_test, samples_for_test = self.reference_data()
        self.assertEqual(sorted(names_from_mdf), sorted(test_channels))
        self.assertEqual(sorted(times_from_mdf), sorted(timestamps_for_test))
        self.assertEqual(sorted(samples_from_mdf), sorted(samples_for_test))
        self.assertEqual(sorted(self.code.read_channels()[1]), sorted(test_channels))

    def read_for_test(self, filepath):

        testmdf = asammdf.mdf.MDF(filepath)
        t_names = []
        t_times = []
        t_samples = []
        iter_test = testmdf.iter_channels(skip_master=True)

        for signal in iter_test:
            t_names.append(signal.name)
            st = signal.timestamps
            ss = signal.samples

            if not isinstance(st, list):
                t_times.append(st.tolist())
            else:
                t_times.append(st)
            if not isinstance(ss, list):
                t_samples.append(ss.tolist())
            else:
                t_samples.append(ss)
        g_names = []
        g_times = []
        g_samples = []
        for index, names in enumerate(t_names):
            if t_samples[index] != [] and t_times[index] != []:
                g_names.append(names)
                g_times.append(t_times[index])
                g_samples.append(t_samples[index])
        return g_names, g_times, g_samples

    def reference_data(self):

        test_channels = ['B', 'C', 'D', 'E', 'A']

        s_timestamps = [numpy.arange(0, 10, 2),
                        numpy.arange(0, 10, 1),
                        numpy.arange(0, 10, 3),
                        numpy.arange(5, 20, 1),
                        numpy.arange(-100, 8, 1)]

        s_samples = [numpy.arange(0, 10, 2),
                     numpy.arange(0, 10, 1),
                     numpy.arange(0, 10, 3),
                     numpy.arange(5, 20, 1),
                     numpy.arange(-100, 8, 1)]

        timestamps_for_test = []
        for time in s_timestamps:
            if not isinstance(time, list):
                t = time.tolist()
                timestamps_for_test.append(t)
            else:
                timestamps_for_test.append(time)

        samples_for_test = []
        for sample in s_samples:
            if not isinstance(sample, list):
                s = sample.tolist()
                samples_for_test.append(s)
            else:
                samples_for_test.append(sample)

        return test_channels, timestamps_for_test, samples_for_test

    def generate_comperable_data_from_read(self, my_mdf):
        names_from_mdf = []
        times_from_mdf = []
        samples_from_mdf = []

        testiter = my_mdf.iter_channels(skip_master=True)
        for signal in testiter:
            names_from_mdf.append(signal.name)
            time = []
            samp = []
            for t in signal.timestamps:
                time.append(int(t))
            for sa in signal.samples:
                samp.append(int(sa))
            times_from_mdf.append(time)
            samples_from_mdf.append(samp)

        return names_from_mdf, times_from_mdf, samples_from_mdf
