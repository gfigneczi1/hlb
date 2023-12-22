#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 16:52:14 2020
@author: jez2bp
"""


import os
import re
import _converter.conversion.pltexpressions


class PltParser:
    """ The PltParser implements the reading of the plt file statements into a
    well structured python list

    Methods:
        -read(path)
        -read_text(path)
        -get_statements(text_lines)
        -parse(line)
        -is_valid_statement(line)
        -is_expression(line)
        -is_valid_expression(line)
        -is_comment(line)
        -get_source_name(line)
        -get_target_name(line)
    """

    def __init__(self):
        # strings which represents comment lines in the plt
        self.comment_signs = ['/', '~']
        # a class which can evaluate the plt expressions
        self.expressions = conversion.pltexpressions.PltExpressions()

    def read(self, path):
        """ The main function of the class which runs the parseing process.
        Only takes one argument which is the path of the plt file.
        Returns with a list of dictionaries.
        """

        text_lines = self.read_text(path)
        statements = self.get_statements(text_lines)
        return statements

    def read_text(self, path):
        """ The function reads a text file with the given path and
        returns a line by line list without the blank lines.
        """

        result = None
        if isinstance(path, str):
            if os.path.isfile(path):
                text_lines = []
                with open(path, 'rb') as file:
                    for rawline in file:
                        text_line = ''
                        try:
                            text_line = rawline.decode('utf-8').replace('\n', '').replace('\r', '')
                        except UnicodeDecodeError:
                            text_line = rawline.decode('latin-1').replace('\n', '').replace('\r', '')
                        if text_line != '':
                            text_lines.append(text_line)
                file.close()
                result = text_lines
            else:
                print("PLT file is not exists: ", path)
                result = False
        else:
            print("Invalid PLT file name: ", path)
        return result

    def get_statements(self, text_lines):
        """ The function takes a list of lines as argument and returns with
        the final parsed data.
        """

        statements = []
        last_statement = None
        last_status = None
        last_type = None
        for rawline in text_lines:
            line = rawline.rstrip()
            # get the neccessary data and the category of the data
            statement, statement_type, status = self.parse(line)
            # if the last_statement is a valid statement we assign the parsed
            # line data to the statements list and set the actual statement to
            # the last_statement for the next iteration
            if statement_type == 'statement':
                if last_statement is not None and last_status:
                    statements.append(last_statement)
                last_statement = statement
                last_type = statement_type
                last_status = status
            # if the parsed line is a valid expression we assign it to the
            # previous line which is acutally the last statement
            elif statement_type == 'expression':
                if last_status and last_type == 'statement':
                    if last_statement is not None:
                        last_statement['exp'] = statement
                else:
                    print("Expression '", statement, "' skipped because the last statement was invalid")
                last_type = statement_type
                last_status = status
        # add the last line to the statements if valid
        if last_statement:
            statements.append(last_statement)
        return statements

    def parse(self, line):
        """ The function handles a single line and decides the category
        of it. Provides the proper data and the result of the line parseing.
        """
        result = (None, 'empty', False)
        if line:
            if self.is_comment(line):
                result = (None, 'comment', False)
            elif self.is_expression(line):
                if self.is_valid_expression(line):
                    result = (line, 'expression', True)
                else:
                    print("Expression is invalid: ", line)
            elif self.is_valid_statement(line):
                statement = {}
                statement["source"] = self.get_source_name(line)
                statement["target"] = self.get_target_name(line)
                result = (statement, 'statement', True)
            else:
                print("Unknown line type: ", line)
                result = (None, 'unknown', False)
        return result

    def is_valid_statement(self, line):
        """ The function returns with True if the text line contains at least
        one string as the start of the line and the string has minimum one
        alphabetic character, otherwise returns with False.
        """

        result = False
        substrings = line.split(" ")
        if substrings:
            sourcename = substrings[0]
            if sourcename:
                for string in sourcename:
                    if string.isalpha():
                        result = True
            else:
                print("Missing source signal name: ", line)
        else:
            print("Invalid line: ", line)
        return result

    def is_expression(self, line):
        """ The function returns with Ture if the line starts with "+ exp"
        otherwise returns with False.
        """
        return line.startswith("+ exp")
    def is_valid_expression(self, line):
        """ The function returns with True if the expression.createFunction
        external function evaluates as valid otherwise returns with False.
        """

        result = True
        # we need to eliminate the '+ ' characters from the line to provide
        # the right format to the createFunction
        function, signalnamelist, expression = self.expressions.createFunction(line[2:])
        if function is None and signalnamelist is None and expression is None:
            result = False
        return result

    def is_comment(self, line):
        """ The function returns with True if the line starts with any comment
        character otherwise it returns with False
        """

        index = 0
        found = False
        while (index < len(self.comment_signs) and not found):
            found = line.startswith(self.comment_signs[index])
            index += 1
        return found

    def get_source_name(self, line):
        """ The function returns with the first space separated string from
        the line
        """

        return line.split(" ")[0]

    def get_target_name(self, line):
        """ The function returns with the value of the "key=value" text
        expression in the line. If the key is not defined it will returns with
        the first space separated value from the line.
        """

        result = self.get_source_name(line)
        regexp = re.compile(r' key+=(\S*)')
        value_of_search = regexp.search(line)
        if value_of_search is not None:
            result = str(value_of_search.group(1))
        return result
