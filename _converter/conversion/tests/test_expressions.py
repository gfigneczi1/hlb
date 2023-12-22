#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 13:16:44 2019

@author: jez2bp
"""

import unittest
import conversion.expression
import numpy
import asammdf

class TestExpressions(unittest.TestCase):

    def setUp(self):

        timestamp_all = numpy.arange(0,180,0.005)
        sample1 = 2*(numpy.arange(0,180,0.005))
        sample2 = 3*(numpy.arange(0,180,0.005))
        sample3 = 4*(numpy.arange(0,180,0.005))
        
        sig1 = asammdf.mdf.Signal(samples=sample1, timestamps=timestamp_all, name='A')
        sig2 = asammdf.mdf.Signal(samples=sample2, timestamps=timestamp_all, name='B')
        sig3 = asammdf.mdf.Signal(samples=sample3, timestamps=timestamp_all, name='C')

        self.mdf = asammdf.mdf.MDF()
        self.mdf.append([sig1])
        self.mdf.append([sig2])
        self.mdf.append([sig3])
        
        self.pltdict = {'A':{'scmin':0,
                             'scmax':1000000,
                             'name':'X',
                             'unit': 'm',
                             'bit': 0,
                             'nobit': 1,
                             'comment': 'In the mdf the name of the signal is X but we renamed to A',
                             'exp': '+ exp=(x/2)/0.005'},
                        'B':{'scmin':0,
                             'scmax':1000000,
                             'name':'Y',
                             'unit': 's',
                             'bit': 0,
                             'nobit': 1,
                             'comment': 'In the mdf the name of the signal is Y but we renamed to B',
                             'exp': '+ exp(x1=A)=(((x/3)/0.005)+x1)/2'},
                        'C':{'scmin':0,
                             'scmax':1000000,
                             'name':'Z',
                             'unit': 'kg',
                             'bit': 0,
                             'nobit': 1,
                             'comment': 'In the mdf the name of the signal is Z but we renamed to C',
                             'exp': '+ exp(x1=A,x2=B)=(((x/4)/0.005)+x1+x2)/3'}}

        self.code=conversion.expression.Expression(self.pltdict, self.mdf)
        

    def test_applicate_expression(self):
        self.assertIsInstance(self.code.applicate_expression(),
                              type(asammdf.mdf.MDF()),
                              'Expression.applicate_expression: checking return value type')

        my_mdf = self.code.applicate_expression()
        reference = list(numpy.arange(0,180000/5,1))
        iteredgroups = my_mdf.iter_channels(skip_master=True)
        samples = []
        for signal in iteredgroups:
            samples.append(signal.samples)

        self.assertEqual(list(samples[0]),
                         list(samples[1]),
                         'Expression.applicate_expression: checking the samples')

        self.assertEqual(list(samples[0]),
                         list(samples[2]),
                         'Expression.applicate_expression: checking the samples')

        for index, item in enumerate(samples):
            for subindex, subitem in enumerate(item):
                samples[index][subindex] = round(subitem,5)

        self.assertEqual(reference,
                         list(samples[0]),
                         'Expression.applicate_expression: checking if the sample is equal with the reference')

    def test_gen_evaldict(self):

        self.code.read_expressions()
        self.assertIsInstance(self.code.gen_evaldict(),
                              dict,
                              'Expression.gen-evaldict: checking return value type')

        reference = {'A': 'exp=(x/2)/0.005',
                     'B': 'exp(x1=A)=(((x/3)/0.005)+x1)/2',
                     'C': 'exp(x1=A,x2=B)=(((x/4)/0.005)+x1+x2)/3'}

        self.assertEqual(reference, self.code.gen_evaldict(),
                         'Expression.gen-evaldict: checking if the evaldict is equal with the reference')

    def test_calc_new_mdf(self):

        self.code.contruct_sourcedict()
        self.code.read_expressions()

        evaldict = {'A': 'exp=(x/2)/0.005',
                    'B': 'exp(x1=A)=(((x/3)/0.005)+x1)/2',
                    'C': 'exp(x1=A,x2=B)=(((x/4)/0.005)+x1+x2)/3'}

        self.assertIsInstance(self.code.calc_new_mdf(evaldict),
                              type(asammdf.mdf.MDF()),
                              'Expression.calc_new_mdf: checking return value type')

        my_mdf = self.code.calc_new_mdf(evaldict)

        reference = list(numpy.arange(0,180000/5,1))
        iteredgroups = my_mdf.iter_channels(skip_master=True)
        samples = []
        for signal in iteredgroups:
            samples.append(signal.samples)

        self.assertEqual(list(samples[0]),
                         list(samples[1]),
                         'Expression.calc_new_mdf: checking the samples')

        self.assertEqual(list(samples[0]),
                         list(samples[2]),
                         'Expression.calc_new_mdf: checking the samples')

        for index, item in enumerate(samples):
            for subindex, subitem in enumerate(item):
                samples[index][subindex] = round(subitem,5)

        self.assertEqual(reference,
                         list(samples[0]),
                         'Expression.calc_new_mdf: checking if the sample is equal with the reference')
