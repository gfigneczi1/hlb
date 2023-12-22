#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 16:29:03 2019

@author: jez2bp
"""
import conversion.reader.read_mdf

import asammdf
import numpy

import unittest

class TestSignalNameHandler(unittest.TestCase):

    def __init__(self, testName):
        super(TestSignalNameHandler, self).__init__(testName)

    def setUp(self):
        self.mdf_object = asammdf.MDF("tests/testMDF2.mf4")
        self.code = conversion.reader.read_mdf.SignalNameHandler(self.mdf_object)

    def test_find(self):
        signalmaps = {}
        refnames = ['A_1', 'B 2', 'C:3', 'D_4:5 6', 'E', 'F', "nincsbenne"]
        names = ['A_1', 'B_2', 'C_3', 'D_4_5_6', 'E', 'F', "nincsbenne"]
        for name in names:
            signalmaps[name] = self.code.find(name)

        names_from_signalmap = []
        test_signalmaps = []
        test_signalmap_items = []
        for signalmapkey in list(signalmaps.keys()):
            test_signalmap = signalmaps[signalmapkey]
            test_signalmaps.append(test_signalmap)
            names_from_signalmap.append(list(test_signalmap.keys())[0])
            test_signalmap_items.append(test_signalmap[list(test_signalmap.keys())[0]])

        self.assertEqual(names_from_signalmap, refnames)

        test_signalmap_types = []
        for test_signalmap in test_signalmaps:
            test_signalmap_types.append(type(test_signalmap))
        type_set = list(set(test_signalmap_types))
        self.assertTrue(len(type_set) == 1)
        self.assertEqual(type_set[0], dict)

        test_signalmap_item_types = []
        for test_signalmap_item in test_signalmap_items:
            test_signalmap_item_types.append(type(test_signalmap_item))
        type_set2 = list(set(test_signalmap_item_types))
        self.assertTrue(len(type_set2) == 1)
        self.assertEqual(type_set2[0], tuple)

        test_locations = {}
        for test_signalmap in test_signalmaps:
            for test_signalmap_key in list(test_signalmap.keys()):
                test_locations[test_signalmap_key] = self.mdf_object.whereis(test_signalmap_key)

        self.assertEqual(test_locations['A_1'], signalmaps['A_1']['A_1'])
        self.assertEqual(test_locations['B 2'], signalmaps['B_2']['B 2'])
        self.assertEqual(test_locations['C:3'], signalmaps['C_3']['C:3'])
        self.assertEqual(test_locations['D_4:5 6'], signalmaps['D_4_5_6']['D_4:5 6'])
        self.assertEqual(test_locations['E'], signalmaps['E']['E'])
        self.assertEqual(test_locations['F'], signalmaps['F']['F'])
        self.assertEqual(test_locations['nincsbenne'], signalmaps['nincsbenne']['nincsbenne'], ())

    def test_isinsecure(self):

        self.assertTrue(self.code.isinsecure('A_1'))
        self.assertTrue(self.code.isinsecure('B_2'))
        self.assertTrue(self.code.isinsecure('C_3'))
        self.assertTrue(self.code.isinsecure('D_4_5_6'))
        self.assertFalse(self.code.isinsecure('E'))
        self.assertFalse(self.code.isinsecure('F'))
        self.assertFalse(self.code.isinsecure('nincsbenne'))

    def test_insecurename(self):

        self.assertEqual(self.code.insecurename('A_1')['A_1'],
                         self.mdf_object.whereis('A_1'))
        self.assertEqual(self.code.insecurename('B_2')['B 2'],
                         self.mdf_object.whereis('B 2'))
        self.assertEqual(self.code.insecurename('C_3')['C:3'],
                         self.mdf_object.whereis('C:3'))
        self.assertEqual(self.code.insecurename('D_4_5_6')['D_4:5 6'],
                         self.mdf_object.whereis('D_4:5 6'))
