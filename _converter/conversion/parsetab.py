
# parsetab.py
# This file is automatically generated. Do not edit.
# pylint: disable=W,C,R
_tabversion = '3.10'

_lr_method = 'LALR'

_lr_signature = 'leftMULDIVleftADDSUBADD DIV EQUAL EXPLEAD FLOAT ID INTEGER LP MUL PI RP SEP SUB VARIABLEDEF X fCOS fEQ fGT fGTE fIOR fLT fLTE fSINformat : EXPLEAD EQUAL exprformat : EXPLEAD LP definitions RP EQUAL exprdefinitions : definition SEP definitionsdefinitions : definitiondefinition : VARIABLEDEF EQUAL IDexpr : LP expr RPexpr : expr operator2 exprexpr : INTEGER\n            | FLOATexpr : VARIABLEDEFexpr : Xexpr : PIoperator2 : ADD\n            | SUB\n            | MUL\n            | DIVexpr : functionfunction : fSIN LP expr RPfunction : fCOS LP expr RPfunction : fEQ LP expr SEP expr RPfunction : fGT LP expr SEP expr RPfunction : fLT LP expr SEP expr RPfunction : fGTE LP expr SEP expr RPfunction : fLTE LP expr SEP expr RPfunction : fIOR LP expr SEP expr RP'
    
_lr_action_items = {'EXPLEAD':([0,],[2,]),'$end':([1,5,7,8,9,10,11,12,41,42,54,55,62,69,70,71,72,73,74,],[0,-1,-8,-9,-10,-11,-12,-17,-7,-6,-18,-19,-2,-20,-21,-22,-23,-24,-25,]),'EQUAL':([2,23,38,],[3,40,51,]),'LP':([2,3,6,13,14,15,16,17,18,19,20,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[4,6,6,30,31,32,33,34,35,36,37,6,-13,-14,-15,-16,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,]),'INTEGER':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[7,7,7,-13,-14,-15,-16,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,]),'FLOAT':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[8,8,8,-13,-14,-15,-16,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,]),'VARIABLEDEF':([3,4,6,24,25,26,27,28,30,31,32,33,34,35,36,37,39,51,56,57,58,59,60,61,],[9,23,9,9,-13,-14,-15,-16,9,9,9,9,9,9,9,9,23,9,9,9,9,9,9,9,]),'X':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[10,10,10,-13,-14,-15,-16,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,]),'PI':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[11,11,11,-13,-14,-15,-16,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,]),'fSIN':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[13,13,13,-13,-14,-15,-16,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,]),'fCOS':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[14,14,14,-13,-14,-15,-16,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,]),'fEQ':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[15,15,15,-13,-14,-15,-16,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,]),'fGT':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[16,16,16,-13,-14,-15,-16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,]),'fLT':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[17,17,17,-13,-14,-15,-16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,]),'fGTE':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[18,18,18,-13,-14,-15,-16,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,]),'fLTE':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[19,19,19,-13,-14,-15,-16,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,]),'fIOR':([3,6,24,25,26,27,28,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[20,20,20,-13,-14,-15,-16,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,]),'ADD':([5,7,8,9,10,11,12,29,41,42,43,44,45,46,47,48,49,50,54,55,62,63,64,65,66,67,68,69,70,71,72,73,74,],[25,-8,-9,-10,-11,-12,-17,25,25,-6,25,25,25,25,25,25,25,25,-18,-19,25,25,25,25,25,25,25,-20,-21,-22,-23,-24,-25,]),'SUB':([5,7,8,9,10,11,12,29,41,42,43,44,45,46,47,48,49,50,54,55,62,63,64,65,66,67,68,69,70,71,72,73,74,],[26,-8,-9,-10,-11,-12,-17,26,26,-6,26,26,26,26,26,26,26,26,-18,-19,26,26,26,26,26,26,26,-20,-21,-22,-23,-24,-25,]),'MUL':([5,7,8,9,10,11,12,29,41,42,43,44,45,46,47,48,49,50,54,55,62,63,64,65,66,67,68,69,70,71,72,73,74,],[27,-8,-9,-10,-11,-12,-17,27,27,-6,27,27,27,27,27,27,27,27,-18,-19,27,27,27,27,27,27,27,-20,-21,-22,-23,-24,-25,]),'DIV':([5,7,8,9,10,11,12,29,41,42,43,44,45,46,47,48,49,50,54,55,62,63,64,65,66,67,68,69,70,71,72,73,74,],[28,-8,-9,-10,-11,-12,-17,28,28,-6,28,28,28,28,28,28,28,28,-18,-19,28,28,28,28,28,28,28,-20,-21,-22,-23,-24,-25,]),'RP':([7,8,9,10,11,12,21,22,29,41,42,43,44,52,53,54,55,63,64,65,66,67,68,69,70,71,72,73,74,],[-8,-9,-10,-11,-12,-17,38,-4,42,-7,-6,54,55,-3,-5,-18,-19,69,70,71,72,73,74,-20,-21,-22,-23,-24,-25,]),'SEP':([7,8,9,10,11,12,22,41,42,45,46,47,48,49,50,53,54,55,69,70,71,72,73,74,],[-8,-9,-10,-11,-12,-17,39,-7,-6,56,57,58,59,60,61,-5,-18,-19,-20,-21,-22,-23,-24,-25,]),'ID':([40,],[53,]),}

_lr_action = {}
for _k, _v in _lr_action_items.items():
   for _x,_y in zip(_v[0],_v[1]):
      if not _x in _lr_action:  _lr_action[_x] = {}
      _lr_action[_x][_k] = _y
del _lr_action_items

_lr_goto_items = {'format':([0,],[1,]),'expr':([3,6,24,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[5,29,41,43,44,45,46,47,48,49,50,62,63,64,65,66,67,68,]),'function':([3,6,24,30,31,32,33,34,35,36,37,51,56,57,58,59,60,61,],[12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,]),'definitions':([4,39,],[21,52,]),'definition':([4,39,],[22,22,]),'operator2':([5,29,41,43,44,45,46,47,48,49,50,62,63,64,65,66,67,68,],[24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,]),}

_lr_goto = {}
for _k, _v in _lr_goto_items.items():
   for _x, _y in zip(_v[0], _v[1]):
       if not _x in _lr_goto: _lr_goto[_x] = {}
       _lr_goto[_x][_k] = _y
del _lr_goto_items
_lr_productions = [
  ("S' -> format","S'",1,None,None,None),
  ('format -> EXPLEAD EQUAL expr','format',3,'p_format1','pltexpressions.py',134),
  ('format -> EXPLEAD LP definitions RP EQUAL expr','format',6,'p_format2','pltexpressions.py',138),
  ('definitions -> definition SEP definitions','definitions',3,'p_definitionsn','pltexpressions.py',142),
  ('definitions -> definition','definitions',1,'p_definitions1','pltexpressions.py',146),
  ('definition -> VARIABLEDEF EQUAL ID','definition',3,'p_definition','pltexpressions.py',150),
  ('expr -> LP expr RP','expr',3,'p_exp_inbracket','pltexpressions.py',160),
  ('expr -> expr operator2 expr','expr',3,'p_exp_2tags','pltexpressions.py',164),
  ('expr -> INTEGER','expr',1,'p_exp_number','pltexpressions.py',168),
  ('expr -> FLOAT','expr',1,'p_exp_number','pltexpressions.py',169),
  ('expr -> VARIABLEDEF','expr',1,'p_exp_signals','pltexpressions.py',173),
  ('expr -> X','expr',1,'p_exp_variables','pltexpressions.py',177),
  ('expr -> PI','expr',1,'p_exp_constant','pltexpressions.py',181),
  ('operator2 -> ADD','operator2',1,'p_operator2','pltexpressions.py',185),
  ('operator2 -> SUB','operator2',1,'p_operator2','pltexpressions.py',186),
  ('operator2 -> MUL','operator2',1,'p_operator2','pltexpressions.py',187),
  ('operator2 -> DIV','operator2',1,'p_operator2','pltexpressions.py',188),
  ('expr -> function','expr',1,'p_exp_function','pltexpressions.py',192),
  ('function -> fSIN LP expr RP','function',4,'p_function_sin','pltexpressions.py',196),
  ('function -> fCOS LP expr RP','function',4,'p_function_cos','pltexpressions.py',200),
  ('function -> fEQ LP expr SEP expr RP','function',6,'p_function_eq','pltexpressions.py',204),
  ('function -> fGT LP expr SEP expr RP','function',6,'p_function_gt','pltexpressions.py',208),
  ('function -> fLT LP expr SEP expr RP','function',6,'p_function_lt','pltexpressions.py',212),
  ('function -> fGTE LP expr SEP expr RP','function',6,'p_function_gte','pltexpressions.py',216),
  ('function -> fLTE LP expr SEP expr RP','function',6,'p_function_lte','pltexpressions.py',220),
  ('function -> fIOR LP expr SEP expr RP','function',6,'p_function_ior','pltexpressions.py',224),
]
