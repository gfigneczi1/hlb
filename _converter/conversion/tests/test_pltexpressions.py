#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 18:04:39 2019

@author: jez2bp
"""
import unittest
import conversion.pltexpressions
import numpy

class TestPLTExpressions(unittest.TestCase):

    def setUp(self):

        self.pltsignalnames = ['A_1',
                               'A_2',
                               'B_3',
                               'C_4',
                               'D_5',
                               'E_6',
                               'B_7',
                               'D_8',
                               'E_9']

        self.sourcedict = {'A_1':numpy.arange(0,10,1, dtype=float),
                           'A_2':numpy.arange(0,20,2),
                           'B_3':numpy.arange(0,30,3),
                           'C_4':numpy.arange(0,40,4),
                           'D_5':numpy.arange(0,50,5),
                           'E_6':numpy.arange(0,60,6),
                           'B_7':numpy.arange(0,70,7),
                           'D_8':numpy.arange(0,80,8),
                           'E_9':numpy.arange(0,90,9)}

        self.code = conversion.pltexpressions.PltExpressions(self.pltsignalnames)

        self.expression1 = 'exp=x/2'
        self.expression2 = 'exp(x1=A_2)=(x1+x)/2'
        self.expression3 = 'exp(x1=B_3,x2=A_2)=(x1+x2+x)/3'
        self.expression4 = 'exp(x1=B_3,x2=A_2)=(x1+x2+1)/(x+1)'

        self.expression_wrong1 = 'exp(x1=A_999)=x+x1/2'
        self.expression_wrong2 = 'exp(x1=A_1)=i_wish_you_can_do_some_magic(x+x1/2)**2'

    def test_createFunction(self):

        self.assertIsInstance(self.code.createFunction(self.expression1),
                              tuple)
        self.assertIsInstance(self.code.createFunction(self.expression2),
                              tuple)
        self.assertIsInstance(self.code.createFunction(self.expression3),
                              tuple)
        self.assertIsInstance(self.code.createFunction(self.expression4),
                              tuple)

        self.assertEqual(str(type(self.code.createFunction(self.expression1)[0])),
                         '<class \''+'function'+'\'>')
        self.assertEqual(str(type(self.code.createFunction(self.expression2)[0])),
                         '<class \''+'function'+'\'>')
        self.assertEqual(str(type(self.code.createFunction(self.expression3)[0])),
                         '<class \''+'function'+'\'>')
        self.assertEqual(str(type(self.code.createFunction(self.expression4)[0])),
                         '<class \''+'function'+'\'>')

        self.assertIsInstance(self.code.createFunction(self.expression1)[1],
                              list)
        self.assertIsInstance(self.code.createFunction(self.expression2)[1],
                              list)
        self.assertIsInstance(self.code.createFunction(self.expression3)[1],
                              list)
        self.assertIsInstance(self.code.createFunction(self.expression4)[1],
                              list)

        self.assertEqual(self.code.createFunction(self.expression1)[2],
                         'x[i]/2')
        self.assertEqual(self.code.createFunction(self.expression2)[2],
                         '(signals[\''+'x1'+'\'][i]+x[i])/2')
        self.assertEqual(self.code.createFunction(self.expression3)[2],
                         '(signals[\''+'x1'+'\'][i]+signals[\''+'x2'+'\'][i]+x[i])/3')
        self.assertEqual(self.code.createFunction(self.expression4)[2],
                         '(signals[\''+'x1'+'\'][i]+signals[\''+'x2'+'\'][i]+1)/(x[i]+1)')

        self.assertIsInstance(self.code.createFunction(self.expression_wrong1),
                              tuple)
        self.assertIsInstance(self.code.createFunction(self.expression_wrong2),
                              tuple)

        self.assertEqual(self.code.createFunction(self.expression_wrong1)[0],
                         None)
        self.assertEqual(self.code.createFunction(self.expression_wrong2)[0],
                         None)

        self.assertEqual(self.code.createFunction(self.expression_wrong1)[1],
                         None)
        self.assertEqual(self.code.createFunction(self.expression_wrong2)[1],
                         None)

        self.assertEqual(self.code.createFunction(self.expression_wrong1)[2],
                         None)
        self.assertEqual(self.code.createFunction(self.expression_wrong2)[2],
                         None)

        signalvalues1 ={}
        signalsamples1 = self.sourcedict['A_2']
        signalsamples_res1 = numpy.arange(0,10,1)
        newfunct, _, _ =self.code.createFunction(self.expression1)

        for index, value in enumerate(signalsamples1):
            signalsamples_res1[index] = newfunct(signalsamples1,
                                                 index,
                                                 signalvalues1)
        self.assertEqual(list(signalsamples1/2),
                         list(signalsamples_res1))

        signalvalues2 ={'x1':self.sourcedict['A_2']}
        signalsamples2 = self.sourcedict['A_1']
        signalsamples_res2 = numpy.arange(0,10,1, dtype=float)
        newfunct, signals, text =self.code.createFunction(self.expression2)

        for index, value in enumerate(signalsamples2):
            signalsamples_res2[index] = newfunct(signalsamples2, index, signalvalues2)
        self.assertEqual([0.0, 1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0, 13.5],
                         list(signalsamples_res2))

        signalvalues3 ={'x1':self.sourcedict['B_3'],
                        'x2':self.sourcedict['A_2']}
        signalsamples3 = self.sourcedict['A_1']
        signalsamples_res3 = numpy.arange(0,10,1, dtype=float)
        newfunct, signals, text =self.code.createFunction(self.expression3)

        for index, value in enumerate(signalsamples3):
            signalsamples_res3[index] = newfunct(signalsamples3, index, signalvalues3)
        self.assertEqual(list(self.sourcedict['A_2']),
                         list(signalsamples_res3))

        signalvalues4 ={'x1':self.sourcedict['B_3'],
                        'x2':self.sourcedict['A_2']}
        signalsamples4 = self.sourcedict['A_1']
        signalsamples_res4 = numpy.arange(0,10,1, dtype=float)
        newfunct, signals, text =self.code.createFunction(self.expression3)

        for index, value in enumerate(signalsamples3):
            signalsamples_res4[index] = newfunct(signalsamples4, index, signalvalues4)
        self.assertEqual(list(self.sourcedict['A_2']),
                         list(signalsamples_res4))
