import pandas as pd
import numpy as np
import copy
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum

class driving_situation_window():
    '''
    classes for classification attribute: upstart, turn, stop, speedup, lane_change
    '''
    def __init__(self, start_index: int, end_index: int, data: pd.DataFrame, classification: DrivingEventEnum=None, trip_name: str=None, index: int=None):
        self.start_index = max(0, start_index)                          #time point where event window starts
        self.end_index = min(end_index, data.shape[0] - 1)              #time point where event window ends
        self.length = self.end_index + 1 - self.start_index                     
        self.data = data.iloc[self.start_index:self.end_index + 1].copy()   #pandas DataFrame that includes all of the data in the event window
        self.classification = classification                            #classification assigned during the event based window making process
        self.driving_event_reason = []
        self.parameters = [0,0,0]                                     #parameters we want to learn at the end of the process
        self.detection_log = ""
        self.trip_name = trip_name
        self.index = index

    def get_time_interval(self) -> list:
        if self.data.shape[0]:
            return [self.data["ts_s_0"].iloc[0], self.data["ts_s_0"].iloc[-1]]
        else:
            return []

    def get_name(self) -> str:
        name = ""
        if self.trip_name is not None:
            name += f"{self.trip_name}"
        if self.classification is not None:
            name += f"_{self.classification.value}"
        if self.index is not None:
            name += f"_{self.index}"
        return name

    def __str__(self):
        if len(self.data) == 0:
            return "driving_situation_window instance empty"
        string = ""
        string = f"driving_situation_window instance {self.get_time_interval()}"
        if self.classification:
            string += f":{self.classification.value}"
        string += f", log:{self.detection_log}"
        return string

    def __eq__(self, other):
        return self.start_index == other.start_index and self.end_index == other.end_index

    def __lt__(self, other):
        if self.start_index != other.start_index:
            return self.start_index < other.start_index
        return self.end_index < other.end_index
