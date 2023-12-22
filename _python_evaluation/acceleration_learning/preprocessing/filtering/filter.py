import os
import sys
import logging
import json
import pandas as pd
from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from utils.convert_unit import convert_m_per_s_to_km_per_h, convert_km_per_h_to_m_per_s

def acc_filter_signal_available(signal_profile: dict) -> bool:
    if signal_profile["ACC_active"] == "unavailable":
        return False
    return True

def detect_acc_active_function(data_frame: pd.DataFrame, index: int, signal_profile: dict) -> bool:
    """Check if the Adaptive Cruise Control is active at the given index in the given data frame.

    Args:
        data_frame (pd.DataFrame): a data frame representing a vehicle trip. It must contain the field "ACC_active".
        index (int): an integer that is in the bounds of the data_frame.

    Returns:
        bool: true if the ACC is active, false elsewhere.
    """
    return 0 < data_frame[signal_profile["ACC_active"]].iloc[index]

def acc_is_activated(data_frame: pd.DataFrame, signal_profile: dict) -> bool:
    """Search the entire data frame for a moment when the Adaptive Cruise Control is activated.

    Args:
        data_frame (pd.DataFrame): a data frame representing a vehicle trip. It must contain the field "ACC_active".

    Returns:
        bool: true if the ACC is on at least once, false elsewhere.
    """
    for i_frame in range(data_frame.shape[0]):
        if detect_acc_active_function(data_frame, i_frame, signal_profile):
            return True
    return False

def filter_acc(windows: [driving_situation_window], signal_profile: dict) -> ([driving_situation_window],[driving_situation_window]):
    """Remove the windows where the Adaptive Cruise Control is activated.

    Args:
        windows (driving_situation_window]): a list of windows with the "data" attribute containing the field "ACC_active".

    Returns:
       [driving_situation_window]: a new list of windows without the windows containing ACC activation
       [driving_situation_window]: a new list of windows containing the ACC activation
    """
    filtered_windows = []
    filtered_out_windows = []
    for window in windows:
        if not acc_is_activated(window.data, signal_profile):
            filtered_windows.append(window)
        else:
            filtered_out_windows.append(window)
            log_dict = {
                "name": window.get_name(),
                "time_interval": window.get_time_interval()
            }
            logging.debug(json.dumps(log_dict))
    return filtered_windows, filtered_out_windows

def vehicle_ahead_filter_signals_available(signal_profile: dict) -> bool:
    if signal_profile["Vehicle_Ahead_Status"] == "unavailable":
        return False
    if signal_profile["Vehicle_Ahead_Following_Time"] == "unavailable":
        return False
    return True

def detect_non_free_driving_function(data_frame: pd.DataFrame, index: int, signal_profile: dict) -> bool:
    """Check if there is a vehicle ahead at the given index in the given data_frame.

    Args:
        data_frame (pd.DataFrame): a data frame representing a vehicle trip. It must contain the fields:
            _ signal_profile["Longitudinal_Velocity"]
            _ signal_profile["Vehicle_Ahead_Status"]
            _ signal_profile["Vehicle_Ahead_Following_Time"]
        index (int): an integer that in the bounds of the data_frame.

    Returns:
        bool: true if a vehicle ahead is detected.
    """
    M_PER_S_TO_KM_PER_H_COEF = 3.6
    SPEED_THRESHOLD_WHERE_ABOVE_VEHICLE_AHEAD_TIME_DETECTION_IS_UNRELIABLE_IN_KM_PER_H = 30
    TIME_GAP_BETWEEN_VEHICLE_AHEAD_TOLERANCE_IN_S = 2.0
    vehicle_ahead_is_close = data_frame[signal_profile["Vehicle_Ahead_Following_Time"]].iloc[index] and data_frame[signal_profile["Vehicle_Ahead_Following_Time"]].iloc[index] < TIME_GAP_BETWEEN_VEHICLE_AHEAD_TOLERANCE_IN_S
    speed_is_low = convert_m_per_s_to_km_per_h(data_frame[signal_profile["Longitudinal_Velocity"]].iloc[index]) < SPEED_THRESHOLD_WHERE_ABOVE_VEHICLE_AHEAD_TIME_DETECTION_IS_UNRELIABLE_IN_KM_PER_H
    vehicle_ahead_detected = 0.0 < data_frame[signal_profile["Vehicle_Ahead_Status"]].iloc[index]
    return (vehicle_ahead_detected and vehicle_ahead_is_close) or (speed_is_low and vehicle_ahead_detected)

def free_driving(data_frame: pd.DataFrame, signal_profile: dict) -> bool:
    """Check if the given data_frame is considered to be a free driving trip.

    Args:
        data_frame (pd.DataFrame): a data frame representing a vehicle trip. It is assumed that the interval of time between two values is 1/30 seconds.
            It must contain the fields:
                _ signal_profile["Time"]
                _ signal_profile["Longitudinal_Velocity"]
                _ signal_profile["Vehicle_Ahead_Status"]
                _ signal_profile["Vehicle_Ahead_Following_Time"]

    Returns:
        bool: true if the trip is a free driving trip
        float: the duration of the non free driving situation
    """
    SAMPLE_RATE_PER_SEC = 30
    FLAT_TIME_TOLERANCE_TO_CONSIDER_THE_VEHICLE_AHEAD_IN_S = 5
    DETECTION_RATIO_TO_CONSIDER_THE_VEHICLE_AHEAD = 0.25
    
    vehicle_ahead_detected_frames = 0
    for i_frame in range(data_frame.shape[0]):
        if detect_non_free_driving_function(data_frame, i_frame, signal_profile):
            vehicle_ahead_detected_frames += 1
    
    vehicle_ahead_detected_duration_in_s = vehicle_ahead_detected_frames / SAMPLE_RATE_PER_SEC
    flat_tolerance_is_exceeded = vehicle_ahead_detected_duration_in_s > FLAT_TIME_TOLERANCE_TO_CONSIDER_THE_VEHICLE_AHEAD_IN_S
    total_duration_in_s = data_frame[signal_profile["Time"]].iloc[-1] - data_frame[signal_profile["Time"]].iloc[0]
    ratio_tolerance_is_exceeded = vehicle_ahead_detected_duration_in_s / total_duration_in_s >= DETECTION_RATIO_TO_CONSIDER_THE_VEHICLE_AHEAD
    if flat_tolerance_is_exceeded or ratio_tolerance_is_exceeded:
        return False, vehicle_ahead_detected_duration_in_s
    else:
        return True, vehicle_ahead_detected_duration_in_s

def filter_non_free_driving(windows: [driving_situation_window], signal_profile: dict) -> ([driving_situation_window],[driving_situation_window]):
    """_summary_

    Args:
        windows (driving_situation_window]): a list of driving situation windows with the "data" attribute having interval of time between two values of 1/30 seconds.
            The "data" attribute must contain the fields:
                _ signal_profile["Vehicle_Ahead_Following_Time"]
                _ signal_profile["Longitudinal_Velocity"]
                _ signal_profile["Vehicle_Ahead_Status"]

    Returns:
       [driving_situation_window]: a new list of windows with only the free driving windows
       [driving_situation_window]: a new list of windows containing the non free driving windows
    """
    filtered_windows = []
    filtered_out_windows = []
    for window in windows:
        is_free_drive, non_free_duration = free_driving(window.data, signal_profile)
        if is_free_drive:
            filtered_windows.append(window)
        else:
            log_dict = {
                "name": window.get_name(),
                "time_interval": window.get_time_interval(),
                "non_free_driving_duration_in_s": non_free_duration
            }
            logging.debug(json.dumps(log_dict))
            filtered_out_windows.append(window)
    return filtered_windows, filtered_out_windows

def extract_filtered_data(data_frame: pd.DataFrame, filter_detection_function, signal_profile) -> [[]]:
    """Extract the index interval where the given filter_detection_function removes data

    Args:
        data_frame (pd.DataFrame): a data frame representing a vehicle trip. It must contain the field that the filter_detection_function uses.
        filter_detection_function (f(the_data_frame, index, signal_profile)->bool): a function that returns true when the filter function would filter at the index in the_data_frame.

    Returns:
        [[int,int]]: a new list of list representing a list of index intervals (2 elements) of the periods when the filter applies.
    """
    filtered_data = []
    filtered_interval = []
    previous_state = filter_detection_function(data_frame, 0, signal_profile)
    if previous_state:
        filtered_interval = [0]
    for i_frame in range(data_frame.shape[0]):
        current_state = filter_detection_function(data_frame, i_frame, signal_profile)
        if not previous_state and current_state:
            filtered_interval = [i_frame]
        if previous_state and not current_state:
            filtered_interval.append(i_frame)
            filtered_data.append(data_frame.iloc[filtered_interval[0]:filtered_interval[1]+1].copy())
        previous_state = current_state
    return filtered_data
