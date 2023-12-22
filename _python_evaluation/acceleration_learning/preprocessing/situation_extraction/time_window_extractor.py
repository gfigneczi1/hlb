import pathlib
import os
import json
import sys
from typing import List
import logging
import numpy as np
from collections import defaultdict
import pandas as pd
from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.preprocessing.situation_extraction import lane_change_cutter, speed_up_cutter, slow_down_cutter
from utils.path import root_path, temp_path

#The extracted time sequence window cannot be shorter than 1s to be used (30 frames per second is frequency)
MINIMUM_WINDOW_TIME_LENGTH_IN_FRAMES = 30
# Stop
EMPIRICAL_LEFT_BOUND_STOP_DISTANCE_FROM_STOP_INDEX = 300
EMPIRICAL_RIGHT_BOUND_STOP_DISTANCE_FROM_STOP_INDEX = 50
# Upstart
EMPIRICAL_LEFT_BOUND_STOP_DISTANCE_FROM_UPSTART_INDEX = 50
EMPIRICAL_RIGHT_BOUND_STOP_DISTANCE_FROM_UPSTART_INDEX = 300
# Turn
EMPIRICAL_LEFT_BOUND_STOP_DISTANCE_FROM_TURN_INDEX = 200
EMPIRICAL_RIGHT_BOUND_STOP_DISTANCE_FROM_TURN_INDEX = 200

class TimeWindowExtractor():
    def __init__(self, trip_name: str):
        self.steeringThreshold = 10
        self.steeringThresholdPassed = False
        self.steeringThresholdPassedCounter = 200
        self.trip_name = trip_name

    '''
    Initial Window making methods 
    '''
    def makeWindows(self, df, signal_profile): 
        '''returns dictionary of lists of driving_situation_window objects

        Keyword arguments:
        df -- pandas DataFrame with all of the data from the inputted file
        '''
        # Event windows
        situation_windows = defaultdict(list)
        stops = []
        upstarts = []
        turns = []
        carFollowingStart = []
        carFollowingEnd = []
        prev = [0, 0]

        # # Situations detections
        for i, (v_x, gamma) in enumerate(zip(df[signal_profile["Longitudinal_Velocity"]],df[signal_profile["YawRate"]])):
            # Stop event
            if v_x == 0 and prev[0] != 0:
                stops.append(i)
            # Upstart event
            elif v_x != 0 and prev[0] == 0:
                upstarts.append(i)
            
            # Turn events
            if gamma > self.steeringThreshold and self.steeringThresholdPassed == False:
                self.steeringThresholdPassed = True
                self.steeringThresholdPassedCounter = 200
                turns.append(i)
            elif gamma < (-1)*self.steeringThreshold and self.steeringThresholdPassed == False:
                self.steeringThresholdPassed = True
                self.steeringThresholdPassedCounter = 200
                turns.append(i)
            if self.steeringThresholdPassed == True:
                self.steeringThresholdPassedCounter -= 1
            if self.steeringThresholdPassedCounter == 0:
                self.steeringThresholdPassed = False
            prev = [v_x, gamma]
        
        situation_windows[DrivingEventEnum.SPEED_UP], _, _ = speed_up_cutter.extract_speed_ups(df, self.trip_name, signal_profile)
        situation_windows[DrivingEventEnum.SLOW_DOWN], _, _ = slow_down_cutter.extract_slow_downs(df, self.trip_name, signal_profile)
        situation_windows[DrivingEventEnum.STOP] = self.makeStopWindows(stops, df)
        situation_windows[DrivingEventEnum.UPSTART] = self.makeUpstartWindows(upstarts, df)
        situation_windows[DrivingEventEnum.TURN] = self.makeTurnWindows(turns, df)
        if self.lane_change_signals_available(signal_profile):
            situation_windows[DrivingEventEnum.LANE_CHANGE] = lane_change_cutter.cut_lane_change(df, self.trip_name, signal_profile)

        named_situation_windows = dict()
        for situations, windows in situation_windows.items():
            named_situation_windows[situations] = self.fill_windows_id_attributes(windows)
        return named_situation_windows

    def extract_window(self, trip_df: pd.DataFrame, situation: DrivingEventEnum) -> [driving_situation_window]:
        """Detect situations and generate the list of the corresponding driving_situation_windows"""
        windows: list = []
        if situation == DrivingEventEnum.SPEED_UP:
            windows, _, _ = speed_up_cutter.extract_speed_ups(trip_df, self.trip_name)
        elif situation == DrivingEventEnum.SLOW_DOWN:
            windows, _, _ = slow_down_cutter.extract_slow_downs(trip_df, self.trip_name)
        elif situation == DrivingEventEnum.LANE_CHANGE:
            windows, _, _ = lane_change_cutter.cut_lane_change(trip_df, self.trip_name)
        else:
            pass
        named_windows = self.fill_windows_id_attributes(windows)
        return named_windows
    
    def fill_windows_id_attributes(self, windows: List[driving_situation_window]) -> [driving_situation_window]:
        for window_id, window in enumerate(windows):
            window.trip_name = self.trip_name
            window.id = window_id
        return windows

    def makeStopWindows(self, stops, df):
        '''returns list of driving_situation_window objects (all with the 'stop' classification)

        Keyword arguments:
        stops -- list including the time points where a stop event was detected
        df -- pandas DataFrame with all of the data from the inputted file
        '''
        stop_windows = []
        for stop in stops:
            stop_window = driving_situation_window(stop - EMPIRICAL_LEFT_BOUND_STOP_DISTANCE_FROM_STOP_INDEX,
                stop + EMPIRICAL_RIGHT_BOUND_STOP_DISTANCE_FROM_STOP_INDEX,
                df,
                DrivingEventEnum.STOP,
            )
            if stop_window.length >= MINIMUM_WINDOW_TIME_LENGTH_IN_FRAMES:
                stop_windows.append(stop_window)
        return stop_windows
    
    def makeUpstartWindows(self, upstarts, df):
        '''returns list of lists of driving_situation_window objects (all with the 'upstart' classification)

        Keyword arguments:
        upstarts -- list including the time points where an upstart event was detected
        df -- pandas DataFrame with all of the data from the inputted file
        '''
        upstart_windows = []
        for upstart in upstarts:
            upstart_window = driving_situation_window(upstart - EMPIRICAL_LEFT_BOUND_STOP_DISTANCE_FROM_UPSTART_INDEX,
                upstart + EMPIRICAL_RIGHT_BOUND_STOP_DISTANCE_FROM_UPSTART_INDEX,
                df,
                DrivingEventEnum.UPSTART,
            )
            upstart_windows.append(upstart_window)
        return upstart_windows
    
    def makeTurnWindows(self, turns, df):
        '''returns list of lists of driving_situation_window objects (all with the 'turn' classification)

        Keyword arguments:
        turns -- list including the time points where a turn event was detected
        df -- pandas DataFrame with all of the data from the inputted file
        '''
        turn_windows = []
        for turn in turns:
            turn_window = driving_situation_window(turn - EMPIRICAL_LEFT_BOUND_STOP_DISTANCE_FROM_TURN_INDEX,
                turn + EMPIRICAL_RIGHT_BOUND_STOP_DISTANCE_FROM_TURN_INDEX,
                df,
                DrivingEventEnum.TURN,
            )
            if turn_window.length >= MINIMUM_WINDOW_TIME_LENGTH_IN_FRAMES:
                turn_windows.append(turn_window)
        return turn_windows

    def lane_change_signals_available(self, signal_profile):
        if signal_profile["Turn_Signals"] == "unavailable":
            return False
        elif signal_profile["Right_Lane_Marking_Lateral_Position"] == "unavailable":
            return False
        elif signal_profile["Left_Lane_Marking_Lateral_Position"] == "unavailable":
            return False
        elif signal_profile["Lateral_Acceleration"] == "unavailable":
            return False
        elif signal_profile["Highway_Status"] == "unavailable":
            return False
        else:
            return True

    '''
    saving methods
    '''
    def make_window_dict(self, window, index, window_index, subject_id):
        '''saves the relevant information about the window object into a json file
        
        keyword arguments:
        window -- the driving_situation_window object whose information we want to save onto a json file
        index -- index to find it among all windows with the same classification
        window_index -- index to find it among all windows
        '''
        info_dict = {
            "trip_name": self.trip_name,
            "subject_id": subject_id,
            "window_index": window_index,
            "start_index": int(window.start_index),
            "end_index": int(window.end_index),
            "length": int(window.length),
            "data": f"{window.classification.value}{index}.csv",
            "classification": window.classification.value,
            "classification_index": index,
            "driving_event_reason": window.driving_event_reason,
            "modeling_parameters": window.parameters
        }
        return info_dict

    def save_window_info(self, window_dict_list: List[dict], window_information_file_path: str):
        ''' saves the information of all the windows (overwrites existing files)

        keyword arguments:
        window_dict_list -- list including all the dictionaries that correspond to all the extracted driving_situation_window objects
        window_information_file_path -- the path to the window information file storing the dictionary
        '''        
        with open(window_information_file_path, 'w', encoding="utf-8") as new_file:
            json.dump(window_dict_list, new_file, indent=4)
    
    def get_saved_windows(self, windows_information_file_path: str, input_time_window_directory_path: str) -> dict:
        '''returns the situation windows extracted from the json file

        keyword arguments:
        windows_information_file_path -- the path to the meta information file of the windows
        input_time_window_directory_path -- the path to the extracted time window output directory (with or without map information)
        '''
        with open(windows_information_file_path, 'r') as open_file:
            read_list = json.load(open_file)
            situation_windows = defaultdict(list)
            for i, window_dict in enumerate(read_list):
                data_path = str(os.path.join(input_time_window_directory_path, window_dict["data"]))
                df = pd.read_csv(data_path)
                if df.shape[0] == 0:
                    continue
                for classification in DrivingEventEnum:
                    if window_dict["classification"] == classification.value:
                        window = driving_situation_window(
                            start_index=0,
                            end_index=df.shape[0]-1,
                            data=df,
                            classification=classification,
                            trip_name=self.trip_name,
                            index=i
                        )
                        window.driving_event_reason = window_dict["driving_event_reason"]
                        window.parameters = window_dict["modeling_parameters"]
                        situation_windows[classification].append(window)
                        break
            return situation_windows

    def save_window_data(self, window, index, output_window_directory_path):
        '''saves window data as well as window information to _temp directory
        
        keyword arguments:
        window -- the driving situation window we want to analyze
        index -- the index of the file sorted by type (e.g. the 2nd detected stop)
        output_window_directory_path -- the directory in the _temp directory where these windows will be saved
        '''
        df = window.data
        if not os.path.exists(output_window_directory_path):
            os.makedirs(output_window_directory_path)
        df.to_csv(os.path.join(output_window_directory_path, window.classification.value + str(index) + ".csv")) #need to change this depending on  the filepath
