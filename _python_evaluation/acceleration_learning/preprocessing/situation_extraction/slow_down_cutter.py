import pandas as pd
import numpy as np
from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.preprocessing.filtering.filter import filter_acc, filter_non_free_driving, acc_filter_signal_available, vehicle_ahead_filter_signals_available
from acceleration_learning.preprocessing.situation_extraction.speed_change_tools import backward_search_for_local_extremum, remove_anomalies, merge_close_windows
from utils.convert_unit import convert_km_per_h_to_m_per_s
from utils.math import value_is_in_tolerance

def cut_slow_down(data: pd.DataFrame, trip_name: str, signal_profile: dict) -> ([driving_situation_window], pd.DataFrame):
    """ Detect and sample the slow down situation windows from the given data frame
    
    Args:
        data (pd.DataFrame): a dataFrame containing the signals of the vehicle during a trip. It is assumed that the interval of time between two values is 1/30 seconds.
            It must contain the signals:
                _ signal_profile["Longitudinal_Velocity"]
                _ signal_profile["Longitudinal_Acceleration"]
                _ signal_profile["Time"]
        trip_name (str): the trip name from which the data come from
    
    Returns:
        [driving_situation_window]: the list of the detected and cut slow downs ordered by time. An empty list is returned if there is no data or speed up found
        pd.DataFrame: a table of the stationary velocity (used for debug purpose). An empty data frame is returned if there is no data
    """
    if not data.shape[0]:
        return [], pd.DataFrame(data={})
    
    NUMBER_OF_SAMPLES_PER_SEC = 30
    VELOCITY_GAP_THRESHOLD_TO_DETECT_SPEEDUP_IN_PERCENTAGE = 0.2
    VELOCITY_GAP_THRESHOLD_TO_DETECT_STABILITY_IN_KM_PER_H = 0.75
    LIFESPAN_TO_BE_STABLE_IN_SEC = 3
    STABILITY_ACCELERATION_TOLERANCE_IN_METER_PER_CUBIC_SEC = 1
    SPEED_DIFFERENCE_THRESHOLD_TO_DETECT_SLOW_PACED_SPEEDUP_IN_KM_PER_H = 10
    
    stationary_state_velocity_sum: float = 0.0
    stationary_state_start_index: int = 0
    is_state_stationary: bool = True
    stationary_velocity: float = 0.0
    
    min_velocity_since_stationary_state_start: float = data[signal_profile["Longitudinal_Velocity"]].iloc[0]
    min_velocity_since_stationary_state_start_index: int = 0
    
    start_speedup_state_acceleration: float = 0
    start_index_of_speed_change_window: int = 0
    end_index_of_speed_change_window: int = 0
    
    slow_downs: list = []
    log = ""
    stationary_velocities = []
    for _ in range(data.shape[0]):
        stationary_velocities.append(0)
    
    for i, v_x in enumerate(data[signal_profile["Longitudinal_Velocity"]]):
        actual_acceleration = data[signal_profile["Longitudinal_Acceleration"]][i]

        if is_state_stationary and v_x > 0.0:
            start_of_speed_change_is_detected = False

            stationary_state_velocity_sum += v_x
            stationary_state_samples_number = i - stationary_state_start_index + 1
            stationary_velocity = stationary_state_velocity_sum / stationary_state_samples_number
            velocity_gap_from_stationary_velocity_in_percent = abs(stationary_velocity - v_x) / max(stationary_velocity, v_x)
            
            # if v_x < min_velocity_since_stationary_state_start:
            #     min_velocity_since_stationary_state_start = v_x
            #     min_velocity_since_stationary_state_start_index = i
            # if data[signal_profile["Longitudinal_Velocity"]].iloc[max(0, i-1)] >= stationary_velocity and v_x < stationary_velocity:
            #     historical_min_velocity_since_stationary_state_start = min_velocity_since_stationary_state_start
            #     historical_min_velocity_since_stationary_state_start_index = min_velocity_since_stationary_state_start_index
            #     min_velocity_since_stationary_state_start = v_x
            #     min_velocity_since_stationary_state_start_index = i
            
            if velocity_gap_from_stationary_velocity_in_percent >= VELOCITY_GAP_THRESHOLD_TO_DETECT_SPEEDUP_IN_PERCENTAGE:
                start_of_speed_change_is_detected = True
                start_index_of_speed_change_window = i
                time = data[signal_profile["Time"]].iloc[start_index_of_speed_change_window]
                log = f"velocity gap from stat velocity [{time}]"
            elif not value_is_in_tolerance(actual_acceleration, STABILITY_ACCELERATION_TOLERANCE_IN_METER_PER_CUBIC_SEC):
                # quick change of acceleration
                start_of_speed_change_is_detected = True
                start_index_of_speed_change_window = i
                time = data[signal_profile["Time"]].iloc[start_index_of_speed_change_window]
                log = f"quick change of acceleration [{time}]"
            elif abs(v_x - min_velocity_since_stationary_state_start) > convert_km_per_h_to_m_per_s(SPEED_DIFFERENCE_THRESHOLD_TO_DETECT_SLOW_PACED_SPEEDUP_IN_KM_PER_H):
                # detect slow change of velocity over time
                start_of_speed_change_is_detected = True
                start_index_of_speed_change_window = min_velocity_since_stationary_state_start_index
                time = data[signal_profile["Time"]].iloc[start_index_of_speed_change_window]
                log = f"slow change of velocity [{time}]"
            if start_of_speed_change_is_detected:
                start_of_speed_change_is_detected = False
                is_state_stationary = False
                start_speedup_state_acceleration = actual_acceleration
        stationary_velocities[i] = stationary_velocity
        if not is_state_stationary:
            end_of_speed_change_is_detected = False
            stationary_velocity_threshold_in_m_per_sec = convert_km_per_h_to_m_per_s(VELOCITY_GAP_THRESHOLD_TO_DETECT_STABILITY_IN_KM_PER_H)
            minimum_duration_to_be_stationary_in_index = LIFESPAN_TO_BE_STABLE_IN_SEC * NUMBER_OF_SAMPLES_PER_SEC
            start_index_of_required_duration_to_be_stable = max(0, i - minimum_duration_to_be_stationary_in_index)
            
            if v_x <= 0.0:
                end_index_of_speed_change_window = i
                end_of_speed_change_is_detected = True
                time = data[signal_profile["Time"]].iloc[i]
                log += f", end:speed is zero [{time}]"
            
            if i - start_index_of_speed_change_window >= minimum_duration_to_be_stationary_in_index and not end_of_speed_change_is_detected:
                start_index_of_required_duration_to_be_stable = max(0, i - minimum_duration_to_be_stationary_in_index)
                max_velocity_during_required_stable_interval = max(data[signal_profile["Longitudinal_Velocity"]][start_index_of_required_duration_to_be_stable:i])
                min_velocity_during_required_stable_interval = min(data[signal_profile["Longitudinal_Velocity"]][start_index_of_required_duration_to_be_stable:i])
                velocity_mean = np.mean(data[signal_profile["Longitudinal_Velocity"]][start_index_of_required_duration_to_be_stable:i])
                max_velocity_is_in_stability_tolerance = value_is_in_tolerance(max_velocity_during_required_stable_interval, stationary_velocity_threshold_in_m_per_sec, velocity_mean)
                min_velocity_is_in_stability_tolerance = value_is_in_tolerance(min_velocity_during_required_stable_interval, stationary_velocity_threshold_in_m_per_sec, velocity_mean)
                # in a fixed period after the current data, are the maximum and minimum value in a stable interval
                if max_velocity_is_in_stability_tolerance and min_velocity_is_in_stability_tolerance:
                    # Add settling
                    end_index_of_speed_change_window = i
                    end_of_speed_change_is_detected = True
                    time = data[signal_profile["Time"]].iloc[end_index_of_speed_change_window]
                    log += f", end:velocity is stable for {LIFESPAN_TO_BE_STABLE_IN_SEC}s [{time}]"

            acceleration_sign_changed = actual_acceleration * start_speedup_state_acceleration < 0
            if acceleration_sign_changed and not value_is_in_tolerance(actual_acceleration, STABILITY_ACCELERATION_TOLERANCE_IN_METER_PER_CUBIC_SEC):
                end_index_of_speed_change_window = i
                end_of_speed_change_is_detected = True
                time = data[signal_profile["Time"]].iloc[end_index_of_speed_change_window]
                log += f", end:acceleration changed greatly [{time}]"

            if end_of_speed_change_is_detected:
                time_interval_in_index = [start_index_of_speed_change_window, end_index_of_speed_change_window]
                if data[signal_profile["Longitudinal_Velocity"]][time_interval_in_index[0]] > data[signal_profile["Longitudinal_Velocity"]][time_interval_in_index[1]]:
                    # slow down case
                    if time_interval_in_index[0] != stationary_state_start_index:
                        less = lambda a, b: a < b
                        local_maximum_index: int = backward_search_for_local_extremum(time_interval_in_index[1], data, signal_profile["Longitudinal_Velocity"], less)
                        if time_interval_in_index[0] < local_maximum_index:
                            time_interval_in_index[1] = local_maximum_index
                        
                        window = driving_situation_window(
                            start_index=time_interval_in_index[0],
                            end_index=time_interval_in_index[1],
                            data=data,
                            classification=DrivingEventEnum.SLOW_DOWN,
                            trip_name=trip_name
                        )
                        window.detection_log = log
                        
                        if window.data.shape[0]:
                            slow_downs.append(window)
                        else:
                            del window
                
                is_state_stationary = True
                stationary_state_start_index = time_interval_in_index[1]
                stationary_state_velocity_sum = v_x
                min_velocity_since_stationary_state_start = v_x
                min_velocity_since_stationary_state_start_index = i
                if start_index_of_required_duration_to_be_stable < time_interval_in_index[1]:
                    stationary_velocity = np.mean(data[signal_profile["Longitudinal_Velocity"]][start_index_of_required_duration_to_be_stable:time_interval_in_index[1]])
                else:
                    stationary_velocity = data[signal_profile["Longitudinal_Velocity"]][1]
    
    stationary_velocity_dict: dict = {signal_profile["Time"]: data[signal_profile["Time"]], "stationary_velocity": stationary_velocities}
    stationary_velocity_df: pd.DataFrame = pd.DataFrame(data=stationary_velocity_dict)
    return slow_downs, stationary_velocity_df

def slow_down_is_abnormal(slow_down: driving_situation_window, signal_profile: dict) -> bool:
    """Detect if a slow down that does not fit into the theoretical speed up tendency of a global decrease of velocity. It is supposed to detect inverted bell-shaped slow down. 

    Args:
        slow_down (driving_situation_window): a driving situation windows of a slow down. Its data attribute must contain the field signal_profile["Longitudinal_Velocity"].

    Returns:
        bool: True if the speed is considered as not normal, False elsewhere.
    """
    start_velocity = slow_down.data[signal_profile["Longitudinal_Velocity"]].iloc[0]
    mean_velocity = np.mean(slow_down.data[signal_profile["Longitudinal_Velocity"]])
    end_velocity = slow_down.data[signal_profile["Longitudinal_Velocity"]].iloc[-1]
    velocity_data_tends_to_increase_from_the_start = start_velocity > mean_velocity
    velocity_data_tends_to_increase_to_the_end = mean_velocity > end_velocity
    return not (velocity_data_tends_to_increase_from_the_start and velocity_data_tends_to_increase_to_the_end)

def filter_slow_downs(slow_downs: [driving_situation_window], signal_profile: dict) -> ([driving_situation_window],{}):
    """Split the slow downs window between the healthy windows and the filtered out windows.

    Args:
        slow_downs (driving_situation_window]): a list of slow down windows. For its data attribute, it is assumed that the interval of time between two values is 1/30 seconds. 
            And it must contain the fields:
                _ signal_profile["Longitudinal_Velocity"]
                _ signal_profile["Longitudinal_Acceleration"]
                _ signal_profile["Time"]
                _ signal_profile["ACC_active"]
                _ signal_profile["Vehicle_Ahead_Following_Time"]
                _ signal_profile["Vehicle_Ahead_Status"]
        
    Returns:
        [driving_situation_window]: the sanitized slow down windows
        dict: a dictionary of the removed windows by the filter
    """
    remaining_slow_downs: list = []
    removed_slow_downs: dict = {}
    if acc_filter_signal_available(signal_profile):
        remaining_slow_downs, removed_slow_downs["filter_acc"] = filter_acc(
            slow_downs,
            signal_profile
        )
    if vehicle_ahead_filter_signals_available(signal_profile):
        remaining_slow_downs, removed_slow_downs["filter_non_free_driving"] = filter_non_free_driving(
            remaining_slow_downs,
            signal_profile
        )
    remaining_slow_downs, removed_slow_downs["anomaly"] = remove_anomalies(
        remaining_slow_downs,
        slow_down_is_abnormal,
        signal_profile
    )
    
    return remaining_slow_downs, removed_slow_downs

def extract_slow_downs(data: pd.DataFrame, trip_name: str, signal_profile: dict) -> ([driving_situation_window], pd.DataFrame, {}):
    """Detect the slow downs situation from the given data, and generate, 
        preprocess and filter the resulting situation windows.
    Args:
        data (pd.DataFrame): a panda DataFrame of a vehicle trip. It is 
            assumed that the interval of time between two values is 1/30 
            seconds.
            It must contain the fields:
                _ signal_profile["Longitudinal_Velocity"]
                _ signal_profile["Longitudinal_Acceleration"]
                _ signal_profile["Time"]
                _ signal_profile["ACC_active"]
                _ signal_profile["Vehicle_Ahead_Following_Time"]
                _ signal_profile["Vehicle_Ahead_Status"]
        trip_name (str): the trip name from which the data comes from
    Returns:
        list: a list of slow downs composed of driving_situation_windows 
            that have been detected and preprocessed. An empty list is 
            returned if there is no data or speed change found for a 
            DrivingEventEnum.
        pd.DataFrame: a table of the stationary velocity (for debug purpose).
            An empty data frame is returned if there is no data
        dict: a dictionary of the speed change windows lists that has been
            removed through the preprocessing process by filtering step, by
            DrivingEventEnum (for debug purpose).
    """
    speed_up_cuts: [driving_situation_window]
    stationary_velocities: pd.DataFrame
    speed_up_cuts, stationary_velocities = cut_slow_down(data, trip_name, signal_profile)

    merged_speed_ups: [driving_situation_window] = merge_close_windows(
        speed_up_cuts,
        data,
        signal_profile
    )

    sanitized_speed_ups: [driving_situation_window]
    removed_speed_ups: dict
    sanitized_speed_ups, removed_speed_ups = filter_slow_downs(merged_speed_ups, signal_profile)

    return sanitized_speed_ups, stationary_velocities, removed_speed_ups
