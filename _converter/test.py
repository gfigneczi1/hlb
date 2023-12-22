#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:56:08 2019

@author: jez2bp
"""
import unittest

import conversion.reader.tests.test_read_MDF_ReadMDF
import conversion.reader.tests.test_read_MDF_SignalNameHandler
import conversion.reader.tests.test_read_d97
import conversion.pltparser.tests.test_parse_PLT
import conversion.reader.tests.test_read_PLT
import conversion.tests.test_pltexpressions
import conversion.tests.test_expressions
import conversion.tests.test_convert_main
import libs.tests.test_scripts
import generate_test_files


class RunAllTest:

    def __init__(self):
        self.alltests = unittest.TestSuite()
        self.force_vals = [True, False]
        self.evalsignals_val_ind = [None, {'test1':[{'name':'A_1'}, {'name':'A_2'}]}, {'test2':[{'name':'A_2'}]}]
        self.inputfilename = ['tests/testMDF.mf4', 'tests/testD97.d97']
        self.GenFile = None

    def run(self):
        self.GenFile = generate_test_files.GenTestFiles('tests')
        self.GenFile.create_test_folder()
        self.GenFile.gen_test_mdf()
        self.GenFile.gen_test_plt()
        self.GenFile.gen_test_mat()
        self.add_tests_for_alltests()
        testrunner = unittest.TextTestRunner()
        testrunner.run(self.alltests)

    def add_tests_for_alltests(self):

        self.alltests.addTest(unittest.makeSuite(libs.tests.test_scripts.TestScripts))
        self.alltests.addTest(unittest.makeSuite(conversion.reader.tests.test_read_d97.TestReadD97))
        self.add_read_mdf_test()
        self.alltests.addTest(unittest.makeSuite(conversion.pltparser.tests.test_parse_PLT.TestParsePLT))
        self.alltests.addTest(unittest.makeSuite(conversion.reader.tests.test_read_PLT.TestReadPLT))
        self.alltests.addTest(unittest.makeSuite(conversion.tests.test_pltexpressions.TestPLTExpressions))
        self.alltests.addTest(unittest.makeSuite(conversion.tests.test_expressions.TestExpressions))
        self.add_convert_main_test()

    def add_read_mdf_test(self):
        self.alltests.addTest(conversion.reader.tests.test_read_MDF_ReadMDF.TestReadMDF('test_read_channels'))
        self.alltests.addTest(conversion.reader.tests.test_read_MDF_SignalNameHandler.TestSignalNameHandler('test_find'))
        self.alltests.addTest(conversion.reader.tests.test_read_MDF_SignalNameHandler.TestSignalNameHandler('test_isinsecure'))
        self.alltests.addTest(conversion.reader.tests.test_read_MDF_SignalNameHandler.TestSignalNameHandler('test_insecurename'))

    def add_convert_main_test(self):
        for file in self.inputfilename:
            for force in self.force_vals:
                    for evalsignals in self.evalsignals_val_ind:
                            self.alltests.addTest(conversion.tests.test_convert_main.TestConvertMain(file, force, evalsignals,
                                                                                                     'test_start_conversion'))
                            self.alltests.addTest(conversion.tests.test_convert_main.TestConvertMain(file, force, evalsignals,
                                                                                                     'test_change_names_and_interpolate'))
                            self.alltests.addTest(conversion.tests.test_convert_main.TestConvertMain(file, force, evalsignals,
                                                                                                     'test_generate_new_signals_and_rename'))
                            self.alltests.addTest(conversion.tests.test_convert_main.TestConvertMain(file, force, evalsignals,
                                                                                                     'test_get_plt_data'))
                            self.alltests.addTest(conversion.tests.test_convert_main.TestConvertMain(file, force, evalsignals,
                                                                                                     'test_get_data'))
                            self.alltests.addTest(conversion.tests.test_convert_main.TestConvertMain(file, force, evalsignals,
                                                                                                     'test_remove_extension'))
                            self.alltests.addTest(conversion.tests.test_convert_main.TestConvertMain(file, force, evalsignals,
                                                                                                     'test_get_dat_mdf_data'))
                            self.alltests.addTest(conversion.tests.test_convert_main.TestConvertMain(file, force, evalsignals,
                                                                                                     'test_get_d97_data'))


RunTests = RunAllTest()
RunTests.run()
