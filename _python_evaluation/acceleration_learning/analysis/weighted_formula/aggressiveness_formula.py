"""Python module to compute a weighted formula that gives a value that represents aggressiveness 
property of a speedup"""
import logging
import numpy as np
import pandas as pd
from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.exception.empty_window_exception import EmptyWindowsException
from utils.math import min_max_normalize

#global variables
## this value was extracted through data analysis (may need to be changed)
### jerk is in m/s^3
LONGITUDINAL_JERK_MAXIMUM_IN_METERS_PER_CUBIC_SECOND = 0.714
LONGITUDINAL_JERK_MINIMUM_IN_METERS_PER_CUBIC_SECOND = -0.792
## These minimum and maximum parameters are taken from the canDB
### this is in %/s
BRAKE_PEDAL_POS_GRADIENT_MAXIMUM_IN_PERCENTAGE_PER_SECOND = 5080
BRAKE_PEDAL_POS_GRADIENT_MINIMUM_IN_PERCENTAGE_PER_SECOND = -5120
### this is in %
ACCELERATION_PEDAL_MAXIMUM_IN_PERCENTAGE = 100
ACCELERATION_PEDAL_MINIMUM_IN_PERCENTAGE = 0

LONGITUDINAL_VELOCITY_MAX_IN_M_PER_S = 127.9375
LONGITUDINAL_VELOCITY_MIN_IN_M_PER_S = -128

ACCELERATION_MAX_IN_M_PER_S_CUBIC = 20.47
ACCELERATION_MIN_IN_M_PER_S_CUBIC = -20.48

def prepare_speedup_formula_parameters(data: pd.DataFrame, signal_profile: dict):
    """Compute the aggressiveness of the data given that it represents a speedup"""
    TEN_MINUTES_IN_S = 60 * 10

    velocity = data[signal_profile["Longitudinal_Velocity"]].to_list()
    acceleration = np.diff(velocity)
    accelerator_pedal_positions = data[signal_profile["Acceleration_Pedal"]].to_list()
    time_in_s = data[signal_profile["Time"]].to_list()

    # Values
    velocity_jerk = np.diff(acceleration)
    norm_max_velocity_jerk = min_max_normalize(
        max(velocity_jerk),
        LONGITUDINAL_JERK_MINIMUM_IN_METERS_PER_CUBIC_SECOND,
        LONGITUDINAL_JERK_MAXIMUM_IN_METERS_PER_CUBIC_SECOND
    )

    norm_max_accelerator_pedal_positions = min_max_normalize(
        max(accelerator_pedal_positions),
        ACCELERATION_PEDAL_MINIMUM_IN_PERCENTAGE,
        ACCELERATION_PEDAL_MAXIMUM_IN_PERCENTAGE
    )

    end_velocity = np.mean(velocity[-data[signal_profile["Longitudinal_Velocity"]].idxmax():])
    velocity_overshoot = max(velocity) - end_velocity
    norm_velocity_overshoot = min_max_normalize(
        velocity_overshoot,
        LONGITUDINAL_VELOCITY_MIN_IN_M_PER_S,
        LONGITUDINAL_VELOCITY_MAX_IN_M_PER_S
    )

    duration_in_s = time_in_s[-1] - time_in_s[0]
    norm_duration = min_max_normalize(duration_in_s, 0, TEN_MINUTES_IN_S)

    parameters = [
        norm_max_velocity_jerk,
        norm_max_accelerator_pedal_positions,
        norm_velocity_overshoot,
        norm_duration
    ]
    return parameters

def evaluate_aggressiveness(
    window: driving_situation_window,
    weights_configuration: dict,
    signal_profile: dict
    ):
    """Compute the aggressiveness of the given driving situation and return it with its parameters
    Raise an EmptyWindowsException if the window's length is less than 3"""
    if window is None or len(window.data) < 3:
        raise EmptyWindowsException
    data = window.data
    dict_weights = weights_configuration[window.classification.value]
    weights = []
    parameters = []

    if window.classification == DrivingEventEnum.SPEED_UP:
        parameters = prepare_speedup_formula_parameters(data, signal_profile)
        weights = [
            dict_weights["velocity_jerk_weight"],
            dict_weights["accelerator_pedal_positions_weight"],
            dict_weights["velocity_overshoot_weight"],
            dict_weights["duration_weight"]
        ]
    
    parameters = np.array(parameters)
    weights = np.array(weights)
    level_of_aggressiveness_dynamics = parameters @ weights

    return level_of_aggressiveness_dynamics, parameters
