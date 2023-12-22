import pandas as pd
import numpy as np
import os
import pathlib

from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
'''
This class assigns reasons for any given acceleration behavior to each window.
The following reasons are considered:
    - stop_sign
    - pedestrian_crossing
    - traffic_light
    - speed_limit_change
    - roundabout
    - highway entrance 
    - highway exit
    - curve in road #not implemented
    - turn 
    - slope #not implemented
'''
MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S = 50 / 3.6
HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S = 80 / 3.6
LIMIT_UNIT = "_miles_per_hour"
YAW_RATE_THRESHOLD_DEG_PER_S = 10

class driving_event_reason_analysis():


    def __init__(self):
        self.stop_sign_present = False
        self.pedestrian_crossing_present = False
        self.traffic_light_present = False
        self.speed_limit_change_present = False
        self.roundabout_present = False
        self.highway_entrance_present = False
        self.highway_exit_present = False
        self.highway_present = False
        self.yaw_rate_threshold_passed = False
    
    def execute(self, window, signal_profile, use_map=False):
        df = window.data
        n = len(df[signal_profile["Longitudinal_Velocity"]])
        
        for i in range(1,n-1):
            if use_map:
                if df["stop_sign_nearby_binary_indicator"][i] == True and not self.stop_sign_present and "stop_sign" not in window.driving_event_reason:
                    self.stop_sign_present = True
                    window.driving_event_reason.append("stop_sign")
                if df["pedestrian_crossing_nearby_binary_indicator"][i] == True and not self.pedestrian_crossing_present and "pedestrian_crossing" not in window.driving_event_reason:
                    self.pedestrian_crossing_present = True
                    window.driving_event_reason.append("pedestrian_crossing")
                if df["traffic_light_nearby_binary_indicator"][i] == True and not self.traffic_light_present and "traffic_light" not in window.driving_event_reason:
                    self.traffic_light_present = True
                    window.driving_event_reason.append("traffic_light")
                if df["roundabouts"][i] == 1 and self.roundabout_present and "roundabout" not in window.driving_event_reason:
                    self.roundabout_present = True
                    window.driving_event_reason.append("roundabout")
                if df["ramps"][i-1] == 1 and df["ramps"][i] == 0:
                    if df["motorway"][i-1] == 0 and df["motorway"][i] == 1 and "highway_entrance" not in window.driving_event_reason:
                        self.highway_entrance_present = True 
                        window.driving_event_reason.append("highway_entrance")
                if df["ramps"][i-1] == 0 and df["ramps"][i] == 1:
                    if df["motorway"][i-1] == 1 and df["motorway"][i] == 0:
                        self.highway_exit_present = True
                        window.driving_event_reason.append("highway_exit")
                try:
                    if df["speed_limits"][i] > 0:
                        speedup_string = "speed_limit_change_from_" + str(df["speed_limits"][i-1]) +"_to_" + str(df["speed_limits"][i] + LIMIT_UNIT)
                        if df["speed_limits"][i-1] != df["speed_limits"][i] and speedup_string not in window.driving_event_reason:
                            self.speed_limit_change_present = True
                            window.driving_event_reason.append(speedup_string)
                except:
                    pass
            if signal_profile["Highway_Status"] != "unavailable":
                #This step only happens if the "Highway_Status" is available in the Signal Profile
                if df[signal_profile["Highway_Status"]].values[i] == 1 and not self.highway_present:
                    self.highway_present = True
                    (window.driving_event_reason.append("highway_driving")) if "highway_driving" not in window.driving_event_reason else None
            if df[signal_profile["YawRate"]].values[i] > YAW_RATE_THRESHOLD_DEG_PER_S and not self.yaw_rate_threshold_passed:
                self.yaw_rate_threshold_passed = True
                (window.driving_event_reason.append("turn")) if "turn" not in window.driving_event_reason else None
        ## add turn detection in here too
        if window.classification == DrivingEventEnum.SPEED_UP:
            long_velocity = np.array(df[signal_profile["Longitudinal_Velocity"]].values)
            min_v = min(long_velocity)
            max_v = max(long_velocity)
            # velocity separation based on min and max values
            category = "speedup_not_categorized_by_velocity"
            if min_v < 0.0:
                category = "backwards_driving"
            elif min_v == 0.0 and max_v < MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S:
                category = "zero_to_low_speed_up"
            elif min_v == 0.0 and max_v > MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S and max_v < HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S:
                category = "zero_to_mid_speed_up"
            elif min_v == 0.0 and max_v > HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S:
                category = "zero_to_high_speed_up"
            elif min_v < MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S and max_v < MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S:
                category = "low_to_low_speed_up"
            elif min_v < MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S and max_v > MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S and max_v < HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S:
                category = "low_to_mid_speed_up"
            elif min_v < MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S and max_v > HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S:
                category = "low_to_high_speed_up"
            elif min_v > MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S and min_v < HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S and max_v > MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S and max_v < HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S:
                category = "mid_to_mid_speed_up"
            elif min_v > MEDIUM_SPEED_VELOCITY_THRESHOLD_M_PER_S and min_v < HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S and max_v > HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S:
                category = "mid_to_high_speed_up"
            elif min_v > HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S and max_v > HIGH_SPEED_VELOCITY_THRESHOLD_M_PER_S:
                category = "high_to_high_speed_up"
            if category not in window.driving_event_reason:
                window.driving_event_reason.append(category)
            if "highway_driving" not in window.driving_event_reason:
                window.driving_event_reason.append("not_highway_driving")
        return window
