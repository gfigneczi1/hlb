"""Module for discrimating driver relatively to others drivers"""
from collections import defaultdict
import numpy as np
import pandas as pd
from acceleration_learning.modeling.model_parameter import CAPT2ModelParameter
from acceleration_learning.analysis.classification.velocity_class import VelocityClass

class Discriminator:
    """ Separate drivers in 2 categories """
    def __init__(self, database: pd.DataFrame):
        self.separators = defaultdict(dict)
        negative_aggressiveness_parameters = [
            CAPT2ModelParameter.D,
            CAPT2ModelParameter.TC,
            CAPT2ModelParameter.TP
        ]
        for velocity_class in list(VelocityClass):
            for model_parameter in list(CAPT2ModelParameter):
                database_where_cur_velocity_class = database[database["velocity_class"] == velocity_class.value]
                if len(database_where_cur_velocity_class) == 0:
                    self.separators[velocity_class][model_parameter] = None
                else:
                    coef = 1
                    if model_parameter in negative_aggressiveness_parameters:
                        coef = -1
                    separation_parameters = {
                        "mean": np.mean(database_where_cur_velocity_class[model_parameter.value]),
                        "std": np.std(database_where_cur_velocity_class[model_parameter.value]),
                        "sign": coef
                    }
                    self.separators[velocity_class][model_parameter] = separation_parameters

    def discriminate(self, value: float, velocity_class: VelocityClass, model_parameter: CAPT2ModelParameter) -> int:
        """Return the category of driving behavior: -1 for calm and 1 for aggressive"""
        if self.separators[velocity_class][model_parameter] is None:
            raise RuntimeError(
                f"there is not separator for the {velocity_class} velocity class with the {model_parameter} model parameter"
            )
        mean = self.separators[velocity_class][model_parameter]["mean"]
        std = self.separators[velocity_class][model_parameter]["std"]
        sign = self.separators[velocity_class][model_parameter]["sign"]
        # 0
        is_on_right_of_distribution = value > mean
        if is_on_right_of_distribution:
            return 1 * sign
        else:
            return -1 * sign
