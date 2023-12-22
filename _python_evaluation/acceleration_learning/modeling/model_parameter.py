"""Module to represent the name of the CA-PT2 model parameters"""
from enum import Enum

class CAPT2ModelParameter(Enum):
    """The Constant Acceleration PT2 model used to model the speed up"""
    D = "D"
    W = "w0"
    TC = "tc"
    A = "alfa"
    TD = "t_delay"
    TP = "Tp"
    OS = "OS"
    TS = "Ts"
