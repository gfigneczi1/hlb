"""Module to profile the user based on model parameters and velocity class"""
from typing import List
import numpy as np
import pandas as pd
from acceleration_learning.modeling.model_parameter import CAPT2ModelParameter
from acceleration_learning.analysis.classification.velocity_class import VelocityClass

def create_profile(driver_id: int, database: pd.DataFrame, discriminator: object) -> List[List[int]]:
    """Return a matrix representing the driver profile horizontally by model parameters and vertically by velocity class"""
    profile = []
    database_where_driver_id = database[database["driver_id"] == driver_id]
    for velocity_class in list(VelocityClass):
        parameters = []
        database_where_driver_id_and_velocity_class = database_where_driver_id[database_where_driver_id["velocity_class"] == velocity_class.value]
        for model_parameter in list(CAPT2ModelParameter):
            select_parameter_from_database_where_driver_id_and_velocity_class = database_where_driver_id_and_velocity_class[model_parameter.value]
            if len(select_parameter_from_database_where_driver_id_and_velocity_class) == 0:
                parameters.append(0)
                continue
            cur_driver_avg = np.mean(select_parameter_from_database_where_driver_id_and_velocity_class)
            parameters.append(discriminator.discriminate(cur_driver_avg, velocity_class, model_parameter))
        profile.append(parameters)
    return profile

def summarize_profile(profile: List[List[int]]) -> pd.DataFrame:
    """Compute an aggressiveness score based on the profile's parameter per velocity class"""
    class_summary = []
    for class_row in profile:
        # D(used)  W0  Tc  a(used) td Tp  Os(used) Ts
        weighted_sum = class_row[0] * 1 +\
            class_row[1] * 0  +\
            class_row[2] * 0 +\
            class_row[3] * 1 +\
            class_row[4] * 0 +\
            class_row[5] * 0 +\
            class_row[6] * 1 +\
            class_row[7] * 0
        weighted_sum /= 3
        class_summary.append(weighted_sum)
    return class_summary