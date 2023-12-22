import sys
import os
import pandas as pd

from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from utils.math import value_is_in_tolerance

def cut_lane_change(data: pd.DataFrame, trip_name: str, signal_profile: dict) -> [driving_situation_window]:
    """ Extract lane change time window from the data

    Args:
        data (pd.DataFrame): a trip data with AVT car signal. It is assumed that the interval of time between two values is 1/30 seconds.
            It must contain the fields:
             _ Longitudinal_Velocity
            _ Right_Lane_Marking_Lateral_Position
            _ Left_Lane_Marking_Lateral_Position
            _ Highway_Status
            _ Lateral_Acceleration
            _ Turn_Signals
        trip_name (str): the trip name from which the data comes from
    Returns:
        [driving_situation_window]: a list of the detected lane changes windows ordered by time
    """
    SAMPLE_RATE_PER_SEC = 30
    MINIMUM_GAP_FROM_LANE_MARKER_TO_BE_INLINE_WITH_THE_ROAD_IN_M = 1
    LAT_ACCELERATION_THRESHOLD_TO_DETECT_LANE_CHANGE = 0.5
    TURN_SIGNAL_OFF = 0
    TURN_SIGNAL_SEARCH_RANGE = SAMPLE_RATE_PER_SEC * 10
    
    lane_is_changing = False
    lane_changes: [driving_situation_window] = []
    # Lane change event
    for i, v_x in enumerate(data[signal_profile["Longitudinal_Velocity"]]):
        right_lane_mark_pos = data[signal_profile["Right_Lane_Marking_Lateral_Position"]].iloc[i]
        left_lane_mark_pos = data[signal_profile["Left_Lane_Marking_Lateral_Position"]].iloc[i]
        if lane_is_changing and left_lane_mark_pos != 0 and right_lane_mark_pos != 0:
            # Lane change has stopped
            lane_is_changing = False
        highway_is_detected = data[signal_profile["Highway_Status"]].iloc[i]
        if v_x > 0 and highway_is_detected and not lane_is_changing and left_lane_mark_pos == 0 and right_lane_mark_pos == 0 and\
            value_is_in_tolerance(data[signal_profile["Lateral_Acceleration"]].iloc[i], LAT_ACCELERATION_THRESHOLD_TO_DETECT_LANE_CHANGE):
            # Lane change has started
            lane_change_start_index = 0
            # Look back when the lane change started, when the turn signal was on
            turn_signal_was_on = False
            backward_search_index = i
            for turn_signal_search_index in range(i, i-TURN_SIGNAL_SEARCH_RANGE, -1):
                if data[signal_profile["Turn_Signals"]].iloc[turn_signal_search_index] != TURN_SIGNAL_OFF:
                    turn_signal_was_on = True
                    backward_search_index = turn_signal_search_index
                    break
            previous_turn_signal = data[signal_profile["Turn_Signals"]].iloc[backward_search_index]
            while backward_search_index >= 0:
                # detect start of lane change if:
                # 1. the turn signal was used: use the turn signal as the start of the lane change
                # 2. (or) the turn signal was not used, but the vehicle was not crossing the line markers
                if turn_signal_was_on:
                    if previous_turn_signal != TURN_SIGNAL_OFF and data[signal_profile["Turn_Signals"]].iloc[backward_search_index] == TURN_SIGNAL_OFF:
                        lane_change_start_index = backward_search_index
                        break
                elif abs(data[signal_profile["Left_Lane_Marking_Lateral_Position"]].iloc[backward_search_index]) > MINIMUM_GAP_FROM_LANE_MARKER_TO_BE_INLINE_WITH_THE_ROAD_IN_M and\
                    abs(data[signal_profile["Right_Lane_Marking_Lateral_Position"]].iloc[backward_search_index]) > MINIMUM_GAP_FROM_LANE_MARKER_TO_BE_INLINE_WITH_THE_ROAD_IN_M:
                        break
                previous_turn_signal = data[signal_profile["Turn_Signals"]][backward_search_index]
                backward_search_index -= 1
            if lane_change_start_index > 0:
                # Look forward when the lane change will stop, when the car is between the two lane markers
                data_length = len(data[signal_profile["Left_Lane_Marking_Lateral_Position"]])
                lane_change_end_index = data_length
                forward_search_index = i
                while forward_search_index < data_length:
                    # detect the end of lane change if:
                    # 1. the both lane markers are detected again and the vehicle is enough centered
                    # 2. (or) at least one lane marker is detected again and teh vehicle is moving laterally
                    if data[signal_profile["Left_Lane_Marking_Lateral_Position"]].iloc[forward_search_index] != 0 and data["RightLaneMrkLatPos_0x347_0"].iloc[forward_search_index] != 0 and\
                    abs(data[signal_profile["Left_Lane_Marking_Lateral_Position"]].iloc[forward_search_index]) > MINIMUM_GAP_FROM_LANE_MARKER_TO_BE_INLINE_WITH_THE_ROAD_IN_M and\
                    abs(data[signal_profile["Right_Lane_Marking_Lateral_Position"]].iloc[forward_search_index]) > MINIMUM_GAP_FROM_LANE_MARKER_TO_BE_INLINE_WITH_THE_ROAD_IN_M:
                        lane_change_end_index = forward_search_index
                        break
                    forward_search_index += 1
                if forward_search_index < data_length:
                    lane_change_window = driving_situation_window(
                        lane_change_start_index,
                        lane_change_end_index,
                        data,
                        DrivingEventEnum.LANE_CHANGE,
                        trip_name=trip_name
                    )
                    lane_changes.append(lane_change_window)
            lane_is_changing = True
    return lane_changes
