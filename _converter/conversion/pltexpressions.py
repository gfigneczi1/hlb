"""
Created on Fri Jan 25 09:56:56 2019
@author: ote2bp
"""
import ply.lex as lex
import ply.yacc as yacc
import math
import re
import logging

class LexerError(Exception): pass
class ParserError(Exception): pass

class PltExpressions:

    lexline = 0

    reserved = {'exp': 'EXPLEAD',
                'SIN': 'fSIN',
                'COS': 'fCOS',
                'EQ': 'fEQ',
                'GT': 'fGT',
                'LT': 'fLT',
                'GTE': 'fGTE',
                'LTE': 'fLTE',
                'IOR': 'fIOR',
                'x' : 'X',
                'PI': 'PI'
}
    tokens = [
        'SEP',
        'LP','RP',
        'EQUAL',
        'INTEGER',
        'FLOAT',
        'ID',
        'VARIABLEDEF',
        'ADD','SUB','MUL','DIV'] + list(reserved.values())
    precedence = (
            ('left','MUL','DIV'),
            ('left', 'ADD','SUB'),
            )

    t_LP = r'\('
    t_RP = r'\)'
    t_EXPLEAD = r'exp'
    t_EQUAL = r'='
    t_SEP = r','
    t_X = r'x'
    t_ADD = r'\+'
    t_SUB = r'\-'
    t_MUL = r'\*'
    t_DIV = r'/'

    def t_PI(self,t):
        r'PI'
        t.value = str(math.pi)
        return t

    #TODO: signs
    def t_FLOAT(self,t):
#        r'[-+]?[0-9]+\.[0-9]+'
        r'[0-9]+\.[0-9]+'
        t.value = t.value
        return t

    #TODO: signs
    def t_INTEGER(self,t):
#        r'[-+]?[0-9]+'
        r'[0-9]+'
        t.value = t.value
        return t

    def t_VARIABLEDEF(self,t):
        r'x[1-9]'
        return t

    def t_ID(self,t):
        r'[A-Za-z_][0-9A-Za-z_]*'
        t.type = self.reserved.get(t.value,'ID')
        return t

    #You should preinit with signalnamelist=["signal1","signal2",...] to check
    #signals are existing in expression
    def __init__(self,signals=None):
        #this makes the parser build
        self.build(debug=0,reflags=re.UNICODE)
        self.signals = signals

    """
        call this for expression evaluation
        returns a lambda function, used signalnamelist, precompiled expression string
        use it on one signal's one value
    """
    def createFunction(self,expression):
        try:
            s1 = self.parse(expression)
            return (eval("lambda x,i,signals :"+s1["expression"]),s1["signals"],s1["expression"])
        except LexerError as err:
            logging.error("{msg} in line {line}".format(msg=err,line=self.error["line"]))
            return (None, None, None)
        except ParserError as err:
            logging.error("{msg} in line {line}".format(msg=err,line=self.error["line"]))
            return (None, None, None)

    #raise a parsererror exception
    def parsererror(self,errordesc,p,tokennum):
        if p != None:
            self.error = {"line":p.lineno(tokennum),"pos":p.lexpos(tokennum)}
        else: self.error = {"line":self.lexline,"pos":0}
        raise ParserError(errordesc)

    #lexical scanner error
    def t_error(self,t):
        self.error = {"line":self.lexline,"pos":0}
        raise LexerError("Illegal character error at '{c}'".format(c=t.value[0]))

    def p_error(self,p):
        if p != None:
            self.parsererror("syntax error at {value}".format(value=p.value),None,1)
        else:
            self.parsererror("syntax error",None,1)

    def build(self,**kwargs):
        self.lexer = lex.lex(module=self, **kwargs)
        self.yaccer = yacc.yacc(module=self)

    #syntacticaL check
    def parse(self,data):
        result = self.yaccer.parse(lexer=self.lexer,input=data)
        return result;

    def p_format1(self,p):
        """format : EXPLEAD EQUAL expr"""
        p[0] = {"expression":p[3],"signals":[]}

    def p_format2(self,p):
        """format : EXPLEAD LP definitions RP EQUAL expr"""
        p[0] = {"expression":p[6],"signals":p[3]}

    def p_definitionsn(self,p):
        "definitions : definition SEP definitions"
        p[0] = [p[1]] + p[3]

    def p_definitions1(self,p):
        "definitions : definition"
        p[0] = [p[1]]

    def p_definition(self,p):
        "definition : VARIABLEDEF EQUAL ID"
        if self.signals != None:
            if p[3] in self.signals:
                p[0] = (p[1],p[3])
            else:
                self.parsererror("invalid signal error:"+p[3],None,1)
        else:
            p[0] = (p[1],p[3])

    def p_exp_inbracket(self,p):
        "expr : LP expr RP"
        p[0] = p[1] + p[2] + p[3]

    def p_exp_2tags(self,p):
        """expr : expr operator2 expr"""
        p[0] = p[1] + p[2] + p[3]

    def p_exp_number(self,p):
        """expr : INTEGER
            | FLOAT"""
        p[0] = p[1]

    def p_exp_signals(self,p):
        "expr : VARIABLEDEF"
        p[0] = "signals['"+p[1]+"'][i]"

    def p_exp_variables(self,p):
        "expr : X"
        p[0] = "x[i]"

    def p_exp_constant(self,p):
        """expr : PI"""
        p[0] = p[1]

    def p_operator2(self,p):
        """operator2 : ADD
            | SUB
            | MUL
            | DIV"""
        p[0] = p[1]

    def p_exp_function(self,p):
        "expr : function"
        p[0] = p[1]

    def p_function_sin(self,p):
        "function : fSIN LP expr RP"
        p[0] = "math.sin("+ p[3] +")"

    def p_function_cos(self,p):
        "function : fCOS LP expr RP"
        p[0] = "math.cos("+ p[3] +")"

    def p_function_eq(self,p):
        "function : fEQ LP expr SEP expr RP"
        p[0] = "(1 if ("+ p[3] + "==" + p[5] +") else 0)"

    def p_function_gt(self,p):
        "function : fGT LP expr SEP expr RP"
        p[0] = "(1 if ("+ p[3] + ">" + p[5] +") else 0)"

    def p_function_lt(self,p):
        "function : fLT LP expr SEP expr RP"
        p[0] = "(1 if ("+ p[3] + "<" + p[5] +") else 0)"

    def p_function_gte(self,p):
        "function : fGTE LP expr SEP expr RP"
        p[0] = "(1 if ("+ p[3] + ">=" + p[5] +") else 0)"

    def p_function_lte(self,p):
        "function : fLTE LP expr SEP expr RP"
        p[0] = "(1 if ("+ p[3] + "<=" + p[5] +") else 0)"

    def p_function_ior(self,p):
        "function : fIOR LP expr SEP expr RP"
        p[0] = "(1 if ("+ p[3] +">0) or ("+ p[5] +">0) else 0)"
