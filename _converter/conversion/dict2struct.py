# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 08:05:52 2021

@author: WIB2BP
"""
##Function desing: list_dict_processor(list input1, list input2)
# Function shall read 2 list: singal names, signal values.
#   \->note: singal names are in the form of nameOfNamespace.nameOfClass.nameOfRunnable.nameOfPort.nameOfSignal.
# Funtion shall get lists as argument.
# Signal name list and signal vaule list shall be aligned.
# The input lists shall be organized into a nested dictionary:
# Dictionary hierarhcy shall be like the following:
# top dict
#  |
#  -- nameOfNamespace dict
#           |
#            -- nameOfClass dict
#                    |
#                     -- nameOfRunnable dict
#                              |
#                              -- nameOfPort dict
#                                     |
#                                     -- nameOfSignal key : signal value
# The valueofsignal = nameOfSignal dict in nameOfPort dict in nameOfRunnable dict in nameOfClass dict in nameOfNamespace.

##Function algo steps
# Iterate trough each member of the singal name list.
#    Parse each signal name member by the specific diveders to get the nameOfNamespace, nameOfClass,...,nameOfSignal
#    Convert the parsed string into nested dictionary.
#    Add the nested dict to the top dict
#    Set the signal value using the list of signal vaule
# Merge nested dicts if the signal_names share common members

##Requirement towards the input:
# The input shall be an ordered list of signals. The ordering key shall be the first common part of the signal name.

class CreateDictStruct:
    
    def __init__(self, signal_name, signal_value):
        self.signal_name = signal_name
        self.signal_value = signal_value
    ##Function definitions
    # The list_dict_processor function takes a list of signal names and values.
    # The function iterates trough on these lists and create a nested dictionary
    # which contains the full signal name as subdicts and the signal value
    def list_dict_processor(self):
        # Declaration of the top dictionary
        top_dict = []
        # Getting length of list and iterating the index
        # same as 'for i in range(len(list))'
        for i in range(len(self.signal_name )):
            # Parsing each member of singal name list
            signal_name_parsed = self.singal_name_parser(self.signal_name[i])
            # Convert each parsed signal name string into a nested dict of the top dict
            nested_dict = self.convert_list_to_dict(signal_name_parsed, self.signal_value[i])
            top_dict.append(nested_dict)
        return top_dict
    #end of def of list_dict_processor function
    
    # The singal_name_parser function take a signal name as a string.
    # The function dived the string by specific devider "." and take every substring into a list.
    # Return the list contain of the substring.
    #TODO: refact signal_name_parser
    def singal_name_parser(self, sigName):
        sgnl_name_parsed = sigName.split(".")
        return sgnl_name_parsed
    #end of def of singal_name_parser function
    
    # The convert_list_to_dict function take a list which
    # contains a parsed signal name and convert it to a nested dict.
    def convert_list_to_dict(self, signal_parsed , signal_val):
        nested_dict = signal_val
        for i in reversed(signal_parsed):
            nested_dict = {i: nested_dict}
        return nested_dict
    #end of def of convert_list_to_dict function
    
    # The merge function is used to merge the nested dicts which contians same signal but different key.
    def merge(self, a, b, path=None):
        "merges b into a"
        if path is None: path = []
        for key in b:
            if key in a:
                if isinstance(a[key], dict) and isinstance(b[key], dict):
                    self.merge(a[key], b[key], path + [str(key)])
                elif a[key] == b[key]:
                    pass # same leaf value
                else:
                    raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
            else:
                a[key] = b[key]
        return a
    # end of def of merge function
    
    # The merge function is used to merge all the nested dicts.
    def merge_all(self, dicts):
        if len(dicts)>=2:
            temp = self.merge(dicts[0], dicts[1])
            for i in range(len(dicts)-2):
                temp = self.merge(temp, dicts[i+2])
        else:
            temp = dicts
        return temp
    # end of def of merge function
    
    # The create_nested_dict function is used to merge all the nested dicts into a top dict.
    def create_nested_dict(self):
        nested_dict = self.list_dict_processor()
        top_nested_dict = self.merge_all(nested_dict)
        return top_nested_dict
    # end of def of create_nested_dict function

