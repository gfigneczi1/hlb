import json
import os
import sys
import pathlib
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from statistics import mean
from enum import Enum

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.preprocessing.situation_extraction.slow_down_cutter import extract_slow_downs
from acceleration_learning.preprocessing.filtering.filter import filter_acc, filter_non_free_driving, extract_filtered_data, detect_acc_active_function, detect_non_free_driving_function

M_PER_S_TO_KM_PER_H_COEF = 3.6
KM_PER_TO_MILE_PER_H_COEF = 0.621371

COEF_FOR_KMH = M_PER_S_TO_KM_PER_H_COEF
COEF_FOR_MPH = M_PER_S_TO_KM_PER_H_COEF * KM_PER_TO_MILE_PER_H_COEF
COEF = COEF_FOR_KMH

def plot_velocity(windows: [driving_situation_window], series_color: str, series_label: str, ax):
    legend_set = False
    for i, window in enumerate(windows):
        print(window)
        data_frame = window.data
        time = data_frame["ts_s_0"]
        velocities = data_frame["LongitudinalVelocity_0x344_20"].apply(lambda x: x*COEF)
        # Plot series
        line, = ax.plot(time, velocities, color=series_color)
        if not legend_set:
            line.set_label(series_label)
            legend_set = True
        # Plot dots
        y = [velocities.iloc[0], velocities.iloc[-1]]
        x = [time.iloc[0], time.iloc[-1]]
        ax.plot(x, y, color=series_color, marker='o', linestyle='None')
        for dot_index in range(len(x)):
            ax.text(x[dot_index], y[dot_index], str(i) + "){:.1f}->{:.1f}".format(x[dot_index], y[dot_index]))

def plot_stationary_velocity(data_frame: pd.DataFrame, ax):
    time = data_frame["ts_s_0"]
    # Plot stationary velocity
    if "stationary_velocity" in data_frame:
        stat_velocities = data_frame["stationary_velocity"].apply(lambda x: x*COEF)
        ax.plot(time, stat_velocities, color="purple", label="Stationary velocity")

# Import the data
expected_timing_file_json: str = f"{pathlib.Path(__file__).resolve().parent.parent.parent}\\tests\\slow_down_moments.json"
expected_timing_file_json = open(expected_timing_file_json, "r")
expected_timings = json.load(expected_timing_file_json)
data_root_folder: str = f"{pathlib.Path(__file__).resolve().parent.parent.parent.parent.parent}\\_temp\\avt_data"

# Compute windows
removed_slowdowns_dict: dict = {}

for trip_key, trip_properties in expected_timings.items():
    filename: str = trip_properties["filename"]
    filepath: str = os.path.join(data_root_folder, filename)
    data_frame: pd.DataFrame = pd.read_csv(filepath)
    cadillac_signal_profile={
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
    slow_downs, stat_velocities, removed_slowdowns_dict = extract_slow_downs(
        data=data_frame,
        trip_name=filename,
        signal_profile=cadillac_signal_profile
    )
    if len(slow_downs):
        ratio = (len(removed_slowdowns_dict["filter_acc"]) + len(removed_slowdowns_dict["filter_non_free_driving"]) + len(removed_slowdowns_dict["anomaly"])) / len(slow_downs)
        print(f"F:{filename}:ratio filter/raws:{ratio}")
    else:
        print(f"F:{filename}: no slow downs")
    
    current_trip_data_frame = data_frame
    
    fig, ax = plt.subplots()
    # Velocity
    window = driving_situation_window(0, len(current_trip_data_frame)-1, current_trip_data_frame)
    plot_velocity([window], "blue", "Vehicle velocity", ax)
    
    # # Expected speed up
    # legend_set = False
    # for expected_speed_up in trip_properties["time"]:
    #     if not legend_set:
    #         ax.axvspan(expected_speed_up[0], expected_speed_up[1], alpha=0.2, color='green', label="Expected speed up periods")
    #         legend_set = True
    #     else:
    #         ax.axvspan(expected_speed_up[0], expected_speed_up[1], alpha=0.2, color='green')
    
    # Slow down
    plot_velocity(slow_downs, "green", "Slow downs", ax)

    # Stationary velocity
    plot_stationary_velocity(stat_velocities, ax)

    # ACC filter
    filtered_data = extract_filtered_data(current_trip_data_frame, detect_acc_active_function, cadillac_signal_profile)
    legend_set = False
    for filter_interval in filtered_data:
        if not legend_set:
            ax.axvspan(filter_interval["ts_s_0"].iloc[0], filter_interval["ts_s_0"].iloc[-1], alpha=0.2, color='red', label='ACC filter periods')
            legend_set = True
        else:
            ax.axvspan(filter_interval["ts_s_0"].iloc[0], filter_interval["ts_s_0"].iloc[-1], alpha=0.2, color='red')
    if "filter_acc" in removed_slowdowns_dict:
        plot_velocity(removed_slowdowns_dict["filter_acc"], "red", "ACC filtered out slow downs", ax)

    # Vehicle ahead filter
    filtered_data = extract_filtered_data(current_trip_data_frame, detect_non_free_driving_function, cadillac_signal_profile)
    legend_set = False
    for filter_interval in filtered_data:
        if not legend_set:
            ax.axvspan(filter_interval["ts_s_0"].iloc[0], filter_interval["ts_s_0"].iloc[-1], alpha=0.2, color='grey', label='Vehicle ahead periods')
            legend_set = True
        else:
            ax.axvspan(filter_interval["ts_s_0"].iloc[0], filter_interval["ts_s_0"].iloc[-1], alpha=0.2, color='grey')
    if "filter_non_free_driving" in removed_slowdowns_dict:
        plot_velocity(removed_slowdowns_dict["filter_non_free_driving"], "grey", "Vehicle ahead filtered out slow downs", ax)

    # Anomalies
    if "anomaly" in removed_slowdowns_dict:
        plot_velocity(removed_slowdowns_dict["anomaly"], "pink", "Anomalies filtered out slow downs", ax)
    
    ax.legend(loc="upper left")
    plt.xlabel("Time (s)")
    if COEF == COEF_FOR_KMH:
        plt.ylabel("Velocity (km/h)")
    elif COEF == COEF_FOR_MPH:
        plt.ylabel("Velocity (mph/h)")
    else:
        plt.ylabel("Velocity ")
    plt.title(f"Velocity graph of the trip named: {filename}")
    plt.show()
