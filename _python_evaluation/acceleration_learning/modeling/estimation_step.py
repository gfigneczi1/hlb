"""Represent the possible step for model parameters estimation"""
from enum import Enum

class EstimationStep(Enum):
    """Model parameter estimation steps"""
    PYTHON_ESTIMATION = "Python_Estimation"
    MATLAB_ESTIMATION_PREPARATION = "Prepare_MATLAB_Estimation"
    MATLAB_COMMANDS_PRINT = "Print_MATLAB_commands"
    MATLAB_RESULTS_SAVE = "Get_MATLAB_results"
