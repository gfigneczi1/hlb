"""Entry file to run a weight formula analysis of speedups"""
import os
import json
import toml
import logging
from acceleration_learning.analysis.weighted_formula.formula_weights_tuning_logic import weighted_formula_tuning_logic
from acceleration_learning.configuration import Configuration, default_configuration
from acceleration_learning.path_configuration import PathConfiguration
from utils.path import root_path

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    configuration = default_configuration()
    weighted_formula_tuning_logic(configuration)
