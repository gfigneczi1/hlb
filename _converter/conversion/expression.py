"""
Created on Thu Jan 24 11:16:01 2019
@author: jez2bp
"""

import re
import logging
import _converter.conversion.pltexpressions as PEX

class Expression:
    """
    The class provides functions to support and handling the data from the .plt
    expression interpreter
    """
    def __init__(self,
                 pltdict,
                 mdf):
        self.mdf = mdf
        self.plt = pltdict
        self.pattern_main = re.compile(r'\+ (.*)')
        self.signal_list = []
        self.expression_list = []
        self.sourcedict = mdf
        self.pltsignalnames = []

    def applicate_expression(self):
        """The main function of the class. Drives the function calling and the
        passing. As result it will returns with an mdf object what contains
        signals with the applied expressions."""
        #self.contruct_sourcedict()
        self.read_expressions()
        evaldict = self.gen_evaldict()
        output_dict = self.calc_new_mdf(evaldict)
        return output_dict

    def contruct_sourcedict(self):
        """The method initialize the source of the samples as a dictionary.
        The key of the dictionray is the name of the signal the value of a key
        is a numoy array with the sample values."""
        iteredgroups = self.mdf.iter_channels(skip_master=True)
        for signal in iteredgroups:
            self.sourcedict[signal.name] = signal.samples

    def read_expressions(self):
        """The method handles the plt data, provides a list with the
        signals and expressions of the plt"""
        sn = []
        exps = []
        for signalname in self.plt:
            self.pltsignalnames.append(signalname)
            if self.plt[signalname]['exp'] is not None:
                sn.append(signalname)
                exps.append(self.plt[signalname]['exp'])
        for index, exp in enumerate(exps):
            if exp != '':
                self.signal_list.append(sn[index])
                self.expression_list.append(exp)

    def gen_evaldict(self):
        """The function make a dictionary for the interpreter method with the
        signal names and expressions"""
        evaluationdict = {}
        for index, item in enumerate(self.expression_list):
            value_s = self.pattern_main.search(item)
            if value_s is not None:
                evaluationdict[self.signal_list[index]] = value_s[1]
        return evaluationdict

    def calc_new_mdf(self, evaluationdict):
        """The function regenerates the new mdf object with the recalculated
        signals with the help of the interpreter method"""
        orderlist = evaluationdict.keys()
        for name in orderlist:
            #iteredgroups = self.mdf.iter_channels(skip_master=True)
            for signalname in self.sourcedict.keys():
                #signalname = signal.name
                if signalname == name:
                    signalsamples = self.sourcedict[signalname]
                    if signalname in evaluationdict.keys():
                        signalvalues = {}
                        expression = evaluationdict[signalname]
                        exp = PEX.PltExpressions(self.pltsignalnames)
                        newfunct, signals, _ = exp.createFunction(expression)
                        try:
                            for sigtokens in signals:
                                if signalname is not None or sigtokens is not None:
                                    signalvalues[sigtokens[0]] = self.sourcedict[sigtokens[1]]
                                else:
                                    logging.error("Unknown expression %s (%s)", str(expression), str(signalname))
                                    return None
                        except Exception:
                            logging.error("Unknown expression %s (%s)", str(expression), str(signalname))
                            return None
                        if newfunct is not None:
                            for index3, _ in enumerate(signalsamples):
                                signalsamples[index3] = newfunct(signalsamples,
                                                                 index3,
                                                                 signalvalues)
                            self.sourcedict[signalname] = signalsamples

        return self.sourcedict
