import json
import os
import sys
import pathlib
from enum import Enum
import argparse
import matplotlib.pyplot as plt
import pandas as pd

from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.preprocessing.situation_extraction.speed_up_cutter import extract_speed_ups, cut_speed_up
from acceleration_learning.preprocessing.situation_extraction.lane_change_cutter import cut_lane_change
from acceleration_learning.preprocessing.filtering.filter import filter_acc, filter_non_free_driving, extract_filtered_data, detect_acc_active_function, detect_non_free_driving_function
from utils.convert_unit import convert_km_per_h_to_mph, convert_m_per_s_to_km_per_h, convert_m_per_s_to_mph
from utils.path import root_path

DEFAULT_DATA_FOLDER: str = f"{root_path()}\\_temp\\avt_data"
DEFAULT_EXPECTED_WINDOWS_FILE_PATH: str = f"{root_path()}\\_python_evaluation\\acceleration_learning\\tests\speed_up_moments.json"

class Metric(Enum):
    M_PER_S = "m_per_s"
    KM_PER_H = "km_per_h"
    MILE_PER_H = "mile_per_h"

def plot_velocity(windows: [driving_situation_window], series_color: str, series_label: str, ax, metric: str, window_bounds_info: bool) -> None:
    legend_set = False
    velocity_converter_coefficient: float = 1.0
    displayed_metric = "m/s"
    if metric == Metric.M_PER_S.value:
        velocity_converter_coefficient = 1.0
        displayed_metric = "m/s"
    elif metric == Metric.KM_PER_H.value:
        velocity_converter_coefficient = convert_m_per_s_to_km_per_h(1.0)
        displayed_metric = "km/h"
    elif metric == Metric.MILE_PER_H.value:
        velocity_converter_coefficient = convert_km_per_h_to_mph(convert_m_per_s_to_km_per_h(1.0))
        displayed_metric = "mph"
        
    for window in windows:
        print(window)
        data_frame = window.data
        time = data_frame["ts_s_0"]
        velocities = data_frame["LongitudinalVelocity_0x344_20"].apply(lambda x: x*velocity_converter_coefficient)
        # Plot velocities series
        line, = ax.plot(time, velocities, color=series_color)
        if not legend_set:
            line.set_label(series_label)
            legend_set = True
        if window_bounds_info:
            # Plot dots for start and end of windows, with annotations
            y = [velocities.iloc[0], velocities.iloc[-1]]
            x = [time.iloc[0], time.iloc[-1]]
            ax.plot(x, y, color=series_color, marker='o', linestyle='None')
            for dot_index in range(len(x)):
                ax.text(x[dot_index] + 2, y[dot_index] + 2, "{:.1f} {}".format(y[dot_index], displayed_metric))
                # ax.text(x[dot_index], y[dot_index], "no." + str(i) + ":{:.1f}->{:.1f}".format(x[dot_index], y[dot_index]))

def display_time_frame(filepath: str, filename: str, show_filtered_period: bool, show_filtered_window: bool, metric: str, expected_windows_dict: dict) -> None:
    # Compute windows
    cadillac_signal_profile = {
        "Time": "ts_s_0",
        "Longitudinal_Velocity": "LongitudinalVelocity_0x344_20",
        "Longitudinal_Acceleration": "ActualAcceleration_0x17D_100",
        "YawRate": "YawRate_0x1E9_20",
        "Highway_Status": "FrwayRoadTypInfo_0x150_100",
        "Acceleration_Pedal": "AcceleratorPedal_0x1A1_25",
        "Longitude": "PsngSysLong_0x32A_100",
        "Latitude": "PsngSysLat_0x32A_100",
        "ACC_active": "ACCactive_0x2CB_0",
        "Vehicle_Ahead_Status": "FOAI_VehicleAhead_0x370_0",
        "Vehicle_Ahead_Following_Time": "VehicleAheadFollowingTimeReq_0x510_100",
        "Vehicle_Ahead_Distance": "VehicleAheadFollowingDistanceReq_0x510_100",
        "Vehicle_Ahead_Velocity": "unavailable",
        "Vehicle_Ahead_Acceleration": "unavailable",
        "Turn_Signals": "TurnSignals_0x140_1000",
        "Lateral_Acceleration": "LateralAcceleration_0x1E9_20",
        "Right_Lane_Marking_Lateral_Position": "RightLaneMrkLatPos_0x347_0",
        "Left_Lane_Marking_Lateral_Position": "LeftLaneMrkLatPos_0x346_0"
    }
    current_trip_data_frame: pd.DataFrame = pd.read_csv(filepath)
    
    speed_ups: list
    removed_speedups_dict: dict
    if show_filtered_window:
        speed_ups, _, removed_speedups_dict = extract_speed_ups(current_trip_data_frame, filename, cadillac_signal_profile)
    else:
        # show only the speed up cut
        speed_ups, _ = cut_speed_up(current_trip_data_frame, filename, cadillac_signal_profile)
    
    if len(speed_ups) == 0:
        print("There is no driving situation windows found")
    
    fig, ax = plt.subplots()
    # Display the whole trip's velocity
    window = driving_situation_window(0, len(current_trip_data_frame)-1, current_trip_data_frame)
    display_info = True
    plot_velocity([window], "blue", "Vehicle velocity", ax, metric, window_bounds_info=False)
    
    # Display the windows
    plot_velocity(speed_ups, "green", "Speed-up windows", ax, metric, display_info)

    # Display the expected time frame of windows 
    for expected_windows_struct in expected_windows_dict.values():
        if expected_windows_struct["filename"] == filename:
            legend_set = False
            for expected_speed_up in expected_windows_struct["time"]:
                if not legend_set:
                    ax.axvspan(expected_speed_up[0], expected_speed_up[1], alpha=0.2, color='green', label="Expected speed-up periods")
                    legend_set = True
                else:
                    ax.axvspan(expected_speed_up[0], expected_speed_up[1], alpha=0.2, color='green')

    # ACC filter
    if show_filtered_period:
        filtered_data = extract_filtered_data(
            current_trip_data_frame,
            detect_acc_active_function,
            cadillac_signal_profile
        )
        legend_set = False
        for filter_interval in filtered_data:
            if not legend_set:
                ax.axvspan(filter_interval["ts_s_0"].iloc[0], filter_interval["ts_s_0"].iloc[-1], alpha=0.2, color='red', label='ACC filter periods')
                legend_set = True
            else:
                ax.axvspan(filter_interval["ts_s_0"].iloc[0], filter_interval["ts_s_0"].iloc[-1], alpha=0.2, color='red')
    if show_filtered_window:
        if "filter_acc" in removed_speedups_dict:
            plot_velocity(removed_speedups_dict["filter_acc"], "red", "ACC filtered out windows", ax, metric, display_info)

    # Vehicle ahead filter
    if show_filtered_period:
        filtered_data = extract_filtered_data(current_trip_data_frame, detect_non_free_driving_function, cadillac_signal_profile)
        legend_set = False
        for filter_interval in filtered_data:
            if not legend_set:
                ax.axvspan(filter_interval["ts_s_0"].iloc[0], filter_interval["ts_s_0"].iloc[-1], alpha=0.2, color='grey', label='Vehicle ahead periods')
                legend_set = True
            else:
                ax.axvspan(filter_interval["ts_s_0"].iloc[0], filter_interval["ts_s_0"].iloc[-1], alpha=0.2, color='grey')
    if show_filtered_window:
        if "filter_non_free_driving" in removed_speedups_dict:
            plot_velocity(removed_speedups_dict["filter_non_free_driving"], "grey", "Vehicle ahead filtered out windows", ax, metric, display_info)

    # Anomalies filtered out
    if show_filtered_window:
        if "anomaly" in removed_speedups_dict:
            plot_velocity(removed_speedups_dict["anomaly"], "pink", "Anomalies filtered out windows", ax, metric, display_info)
    
    ax.legend(loc="upper left")
    plt.xlabel("Time (s)")
    plt.ylabel(f"Velocity ({metric})")
    plt.title(f"Graph of the longitudinal velocity from the trip {filename}")
    plt.grid()
    plt.show()

def plot_marker_distance(windows: [driving_situation_window], right_series_color: str, left_series_color: str, series_label: str, ax, window_bounds_info: bool) -> None:
    legend_set = False
    for i, window in enumerate(windows):
        print(window)
        data_frame = window.data
        time = data_frame["ts_s_0"]
        min_left_marker = min(data_frame["LeftLaneMrkLatPos_0x346_0"])
        left_marker = data_frame["LeftLaneMrkLatPos_0x346_0"].apply(lambda x: x+min_left_marker)
        
        max_right_marker = max(data_frame["RightLaneMrkLatPos_0x347_0"])
        right_marker = data_frame["RightLaneMrkLatPos_0x347_0"].apply(lambda x: x-max_right_marker)
        # Left lane
        left_line, = ax.plot(time, left_marker, color=left_series_color)
        if not legend_set:
            left_line.set_label("Left marking position of " + series_label)
        if window_bounds_info:
            # Plot dots for start and end of windows, with annotations
            y = [left_marker.iloc[0], left_marker.iloc[-1]]
            x = [time.iloc[0], time.iloc[-1]]
            ax.plot(x, y, color=left_series_color, marker='o', linestyle='None')
            for dot_index in range(len(x)):
                ax.text(x[dot_index], y[dot_index], "{:.1f} m".format(y[dot_index]))
        
        # Right lane
        right_line, = ax.plot(time, right_marker, color=right_series_color)
        if not legend_set:
            right_line.set_label("Right marking position of " + series_label)
            legend_set = True
        if window_bounds_info:
            # Plot dots for start and end of windows, with annotations
            y = [right_marker.iloc[0], right_marker.iloc[-1]]
            x = [time.iloc[0], time.iloc[-1]]
            ax.plot(x, y, color=right_series_color, marker='o', linestyle='None')
            for dot_index in range(len(x)):
                ax.text(x[dot_index], y[dot_index], "{:.1f} m".format(y[dot_index]))

def display_lane_change_frame(filepath: str, filename: str):
       # Compute windows
    lane_changes: list
    current_trip_data_frame: pd.DataFrame = pd.read_csv(filepath)
    cadillac_signal_profile = {
        "Time": "ts_s_0",
        "Longitudinal_Velocity": "LongitudinalVelocity_0x344_20",
        "Longitudinal_Acceleration": "ActualAcceleration_0x17D_100",
        "YawRate": "YawRate_0x1E9_20",
        "Highway_Status": "FrwayRoadTypInfo_0x150_100",
        "Acceleration_Pedal": "AcceleratorPedal_0x1A1_25",
        "Longitude": "PsngSysLong_0x32A_100",
        "Latitude": "PsngSysLat_0x32A_100",
        "ACC_active": "ACCactive_0x2CB_0",
        "Vehicle_Ahead_Status": "FOAI_VehicleAhead_0x370_0",
        "Vehicle_Ahead_Following_Time": "VehicleAheadFollowingTimeReq_0x510_100",
        "Vehicle_Ahead_Distance": "VehicleAheadFollowingDistanceReq_0x510_100",
        "Vehicle_Ahead_Velocity": "unavailable",
        "Vehicle_Ahead_Acceleration": "unavailable",
        "Turn_Signals": "TurnSignals_0x140_1000",
        "Lateral_Acceleration": "LateralAcceleration_0x1E9_20",
        "Right_Lane_Marking_Lateral_Position": "RightLaneMrkLatPos_0x347_0",
        "Left_Lane_Marking_Lateral_Position": "LeftLaneMrkLatPos_0x346_0"
    }
    lane_changes = cut_lane_change(current_trip_data_frame, filename, cadillac_signal_profile)
    
    if len(lane_changes) == 0:
        print("There is no driving situation windows found")
    
    fig, ax = plt.subplots()
    # Display the whole trip's velocity
    window = driving_situation_window(0, len(current_trip_data_frame)-1, current_trip_data_frame)
    plot_marker_distance([window], "blue", "orange", "entire trip", ax, False)
    
    # Display the windows
    plot_marker_distance(lane_changes, "green", "green", "Lane change windows", ax, True)

    # # Display the expected time frame of windows 
    # for expected_windows_struct in expected_windows_dict.values():
    #     if expected_windows_struct["filename"] == filename:
    #         legend_set = False
    #         for expected_speed_up in expected_windows_struct["time"]:
    #             if not legend_set:
    #                 ax.axvspan(expected_speed_up[0], expected_speed_up[1], alpha=0.2, color='green', label="Expected windows periods")
    #                 legend_set = True
    #             else:
    #                 ax.axvspan(expected_speed_up[0], expected_speed_up[1], alpha=0.2, color='green')
    
    ax.legend(loc="upper left")
    plt.xlabel("Time (s)")
    plt.ylabel(f"Distance (m)")
    plt.title(f"Graph of the marking positions from the trip named {filename}")
    plt.show() 

def setup_argument_parser() -> argparse.ArgumentParser:
    parser: ArgumentParser = argparse.ArgumentParser(
        description="Display time frame of the speed up window"
    )
    parser.add_argument(
        "--directory",
        dest="data_directory",
        type=str,
        nargs=1,
        default=DEFAULT_DATA_FOLDER,
        help="the path to the trip directory"
    )
    parser.add_argument(
        "--hide_filter_period",
        dest="hide_filtered_period",
        action='store_true',
        default=False,
        help="hide the time period where the filter(s) have effect"
    )
    parser.add_argument(
        "--hide_filtered_window",
        dest="hide_filtered_window",
        action='store_true',
        default=False,
        help="hide the driving window removed by the filter(s)"
    )
    parser.add_argument(
        "--metric",
        dest="velocity_metric",
        type=str,
        nargs=1,
        choices=[Metric.M_PER_S.value, Metric.KM_PER_H.value, Metric.MILE_PER_H.value],
        default=Metric.KM_PER_H.value,
        help="set the metric of the displayed velocity"
    )
    parser.add_argument(
        "--expected_window",
        dest="expected_windows_filepath",
        type=str,
        nargs=1,
        default=DEFAULT_EXPECTED_WINDOWS_FILE_PATH,
        help="the expected timing frame for the windows"
    )
    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument(
        "-f",
        "--filename",
        required=True,
        dest="filename",
        type=str,
        nargs=1,
        help="the name of the trip file"
    )
    required_arguments.add_argument(
        "-s",
        "--scenario",
        required=True,
        dest="scenario",
        type=str,
        nargs=1,
        choices=[DrivingEventEnum.SPEED_UP.value, DrivingEventEnum.LANE_CHANGE.value],
        default=DrivingEventEnum.SPEED_UP.value,
        help="the type of driving situation"
    )
    return parser

if __name__ == "__main__":
    parser: argparse.ArgumentParser = setup_argument_parser()
    args = parser.parse_args()
    expected_timing_file_json = open(args.expected_windows_filepath, "r")
    expected_timings = json.load(expected_timing_file_json)
    filepath: str = os.path.join(args.data_directory[0], args.filename[0])
    if args.scenario[0] == DrivingEventEnum.SPEED_UP.value:
        display_time_frame(
            filepath,
            args.filename[0],
            not args.hide_filtered_period,
            not args.hide_filtered_window,
            args.velocity_metric,
            expected_timings
        )
    elif args.scenario[0] == DrivingEventEnum.LANE_CHANGE.value:
        display_lane_change_frame(
            filepath,
            args.filename[0]
        )
