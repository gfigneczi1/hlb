import copy
import logging
import json
import pandas as pd
from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum

def backward_search_for_local_extremum(start_index: int, data_frame: pd.DataFrame, signal_name: str, comparison_function_for_stop_condition) -> int:
    """Search for the local extremum (minimum or maximum) in the given data frame for the given signal
    
    Args:
        start_index (int): an integer where to start the backward search. It must be in the bound of the given data frame
        data_frame (pd.DataFrame): a data frame that includes the signal_name field in its table.
        signal_name (str): a string of a field in the data_frame
        comparison_function_for_stop_condition (f(a: int, b: int) -> bool): a pointer to a function that compare two values to define how to search for a extremum. The first 'a' parameter is the value with 
            the smaller successive index in the search than the 'b' parameter
    
    Returns:
        int: the index of the local extremum or 0 if there is extremum found
    """
    local_extremum_search_index = start_index
    while local_extremum_search_index - 1 > 0 and not comparison_function_for_stop_condition(data_frame[signal_name][local_extremum_search_index - 1], data_frame[signal_name][local_extremum_search_index]):
        if data_frame[signal_name][local_extremum_search_index] <= 0:
            break
        local_extremum_search_index -= 1
    return local_extremum_search_index

def merge_close_windows(windows: [driving_situation_window], trip_data: pd.DataFrame, signal_profile: dict) -> [driving_situation_window]:
    """Concatenate the windows that are close to each other in time

    Args:
        windows ({driving_situation_window]): a list of speed up situation
            windows, it is assumed that the windows do not overlap in time
        trip_data (pd.DataFrame): the pandas DataFrame of the whole trip 
            where the windows are extracted from. It must includes the data
            of the given windows. 
            It must contain the signal: signal_profile["Time"]

    Returns:
        [driving_situation_window]: a new list of windows with different 
            windows that are the merges of windows
    """
    SHORT_ENOUGH_DURATION_TO_MERGE_WINDOWS_IN_S = 2
    new_windows: [driving_situation_window] = []
    logger = logging.getLogger()
    if len(windows):
        new_windows.append(copy.deepcopy(windows[0]))
        i = 1
        while i < len(windows):
            windows_time_gap_is_short_enough_to_merge = windows[i].data[signal_profile["Time"]].iloc[0] - new_windows[-1].data[signal_profile["Time"]].iloc[-1] <= SHORT_ENOUGH_DURATION_TO_MERGE_WINDOWS_IN_S
            if windows_time_gap_is_short_enough_to_merge:
                last_window = new_windows.pop()
                merged_window = driving_situation_window(
                    start_index=last_window.start_index,
                    end_index=windows[i].end_index,
                    data=trip_data,
                    classification=windows[i].classification,
                    trip_name=windows[i].trip_name
                )
                merged_window.detection_log = f"{last_window.detection_log}|{windows[i].detection_log}"
                new_windows.append(merged_window)
                
                log_dict = {
                    "name": windows[i].get_name(),
                    "merged_windows" : [
                        {
                            "time_interval": last_window.get_time_interval(),
                        },
                        {
                            "time_interval": windows[i].get_time_interval(),
                        }
                    ]
                }
                logger.debug(json.dumps(log_dict))
            else:
                new_windows.append(copy.deepcopy(windows[i]))
            i += 1
    return new_windows

def remove_anomalies(speed_changes: [driving_situation_window], anomaly_detection_function, signal_profile: dict) -> ([driving_situation_window],[driving_situation_window]):
    """Remove the speed changes that are considered as abnormal by the given anomaly_detection_function function

    Args:
        speed_changes (driving_situation_window]): a list of situation windows. The data attribute must contain the fields required by the anomaly_detection_function

    Returns:
        [driving_situation_window]: a new list of speed changes without anomaly
        [driving_situation_window]: a new list of speed changes of the anomalies found
    """
    healthy_speed_changes: [driving_situation_window] = []
    anomalies: [driving_situation_window] = []
    for speed_change in speed_changes:
        if speed_change.data.shape[0]:
            if anomaly_detection_function(speed_change, signal_profile):
                anomalies.append(copy.deepcopy(speed_change))
            else:
                healthy_speed_changes.append(copy.deepcopy(speed_change))
    return healthy_speed_changes, anomalies
