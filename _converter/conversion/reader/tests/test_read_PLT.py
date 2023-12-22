#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 16:50:41 2019

@author: jez2bp
"""
import conversion.reader.read_plt

import unittest

class TestReadPLT(unittest.TestCase):

    def setUp(self):

        self.pltpath = 'tests/testPLT.PLT'
        self.code = conversion.reader.read_plt.ReadPLT(self.pltpath)

    def test_extract_PLT_data(self):

        self.assertIsInstance(self.code.extract_plt_data(),
                              tuple)

        self.assertIsInstance(self.code.extract_plt_data()[0],
                              list)

        self.assertIsInstance(self.code.extract_plt_data()[1],
                              dict)

        test_chan = ['A', 'A', 'B', 'C', 'D', 'E', 'B', 'D' ,'E']

        self.assertEqual(test_chan, self.code.extract_plt_data()[0])
        self.assertEqual(self.dictionary_data_for_testing(), self.code.extract_plt_data()[1])

    def test_readplt(self):

         self.assertIsInstance(self.code.readplt(),
                              dict)

         self.assertEqual(self.code.readplt(), self.dictionary_data_for_testing())

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