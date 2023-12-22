"""Python module analyzing speed up velocity curvature"""
import os
import logging
import numpy as np
import pandas as pd
from utils.math import min_max_normalize
from acceleration_learning.subject_identifier import SubjectIdentifier
from acceleration_learning.path_configuration import PathConfiguration
from acceleration_learning.trip_file import TripFile
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.preprocessing.situation_extraction.time_window_extractor import TimeWindowExtractor

def curvature(speed_up_df: pd.DataFrame, signal_profile: dict) -> float:
    """ Return the area difference between the velocity in the given speed up 
    data frame and a theoretical straight line from the start velocity to the 
    max velocity of the speed up data frame. A positive value means the curve
    is above the straight line, a negative value means the curve is below
    the line."""
    if speed_up_df.empty:
        raise ValueError("Cannot compute curvature if the data frame is empty")
    peak_instant_index = speed_up_df[signal_profile["Longitudinal_Velocity"]].idxmax()
    peak_velocity = speed_up_df[signal_profile["Longitudinal_Velocity"]].iloc[peak_instant_index]
    start_velocity = speed_up_df[signal_profile["Longitudinal_Velocity"]].iloc[0]
    n_samples = peak_instant_index + 1
    straight_line_from_start_to_peak = np.linspace(
        start=start_velocity,
        stop=peak_velocity,
        num=n_samples
    )
    velocity_sequence_from_start_to_peak = speed_up_df[signal_profile["Longitudinal_Velocity"]][0:n_samples]
    diff_sum = 0
    max_velocity = peak_velocity
    min_velocity = start_velocity
    for straight_line_velocity, real_velocity in zip(straight_line_from_start_to_peak, velocity_sequence_from_start_to_peak):
        diff_sum += min_max_normalize(real_velocity, min_velocity, max_velocity) - min_max_normalize(straight_line_velocity, min_velocity, max_velocity)
    area_difference = diff_sum / n_samples
    return area_difference

def menger_curvature(x: list, y: list, z: list) -> float:
    """Return the curvature of a circle having x, y and z as included points"""
    xy_dist = np.sqrt((x[0] - y[0])**2 + (x[1] - y[1])**2)
    yz_dist = np.sqrt((z[0] - y[0])**2 + (z[1] - y[1])**2)
    zx_dist = np.sqrt((z[0] - x[0])**2 + (z[1] - x[1])**2)
    semi_perimeter = (xy_dist + yz_dist + zx_dist) / 2
    triangle_area = np.sqrt(semi_perimeter * (semi_perimeter - xy_dist) * (semi_perimeter - yz_dist) * (semi_perimeter - zx_dist))   # Heron's formula
    if xy_dist == 0 or yz_dist == 0 or zx_dist == 0:
        return 0
    return (4 * triangle_area) / (abs(xy_dist) * abs(yz_dist) * abs(zx_dist))

def curvature_with_circle(data: pd.DataFrame, signal_profile: dict) -> float:
    """Compute the speed up velocity curvature based on the circle formula (Heron)"""
    if data.empty:
        raise ValueError("Cannot compute Heron curvature if the data frame is empty")
    peak_instant_index = data[signal_profile["Longitudinal_Velocity"]].idxmax()
    peak_velocity = data[signal_profile["Longitudinal_Velocity"]].iloc[peak_instant_index]
    start_velocity = data[signal_profile["Longitudinal_Velocity"]].iloc[0]

    max_velocity = peak_velocity
    min_velocity = start_velocity

    curve_start_point = [0, min_max_normalize(start_velocity, min_velocity, max_velocity)]
    curve_middle_index = int(peak_instant_index / 2)
    curve_middle_velocity = data[signal_profile["Longitudinal_Velocity"]].iloc[curve_middle_index]
    curve_middle_point = [
        curve_middle_index,
        min_max_normalize(curve_middle_velocity, min_velocity, max_velocity)
    ]
    curve_end_point = [
        peak_instant_index,
        min_max_normalize(peak_velocity, min_velocity, max_velocity)
    ]
    circle_curvature = menger_curvature(curve_start_point, curve_middle_point, curve_end_point)

    slope = (curve_end_point[1] - curve_start_point[1]) / (curve_end_point[0] - curve_start_point[0])
    middle_velocity_is_below_the_slope = curve_middle_velocity < curve_middle_index * slope
    if middle_velocity_is_below_the_slope:
        circle_curvature = circle_curvature * -1
    return circle_curvature

def compute_speedup_curvature(
    signal_profile: dict,
    path_configuration: PathConfiguration
    ) -> pd.DataFrame:
    """Run the curvature computation using the stored speed ups data and meta information
        Speedups that are longer than 10 seconds are skipped
        Return a table matching speed up, driver and curvature values
    """
    MAXIMUM_SPEEDUP_DURATION_IN_SEC = 10
    speed_up_directory_path: str = path_configuration.preprocessing_extraction_output_directory_path()
    subjects = SubjectIdentifier()
    curvature_table_rows = list()
    n_trip_directory: int = len(os.listdir(speed_up_directory_path))
    for i_trip, trip_directory_name in enumerate(os.listdir(speed_up_directory_path)):
        trip_file = TripFile(os.path.join(path_configuration.preprocessing_input_directory(), trip_directory_name))
        trip_label = '_'.join(trip_directory_name.split("_")[:3])
        extractor = TimeWindowExtractor(trip_name=trip_label)
        situations_windows: dict = extractor.get_saved_windows(
            path_configuration.preprocessing_extraction_information_output_file_path(trip_file),
            os.path.join(speed_up_directory_path, trip_directory_name)
        )
        if not situations_windows.get(DrivingEventEnum.SPEED_UP):
            continue
        driver_id: str = subjects.get_driver_id(trip_label)
        logging.info("Processing the trip: %s (%i/%i)", trip_label, i_trip+1, n_trip_directory)
        for speedup_window in situations_windows[DrivingEventEnum.SPEED_UP]:
            speedup_data_frame = speedup_window.data
            speed_up_duration_in_s = speedup_data_frame[signal_profile['Time']].iloc[-1] - speedup_data_frame[signal_profile['Time']].iloc[0]
            if speed_up_duration_in_s > MAXIMUM_SPEEDUP_DURATION_IN_SEC:
                continue
            speed_up_curvature = curvature(speedup_data_frame, signal_profile)
            speed_up_fullname = trip_label + '-' + str(speedup_window.classification.value) + '-' + str(speedup_window.index)
            curvature_table_rows.append({
                "driver_id": driver_id,
                "speedup_name": speed_up_fullname,
                "curvature": speed_up_curvature
            })
    return pd.DataFrame(curvature_table_rows)

def summary(curvature_df: pd.DataFrame) -> pd.DataFrame:
    """Return table that statistically describe the curvature property from the
    given DataFrame by driver"""
    curvature_summary_rows = []
    for driver_id in pd.unique(curvature_df["driver_id"]):
        curvature_from_driver_id = curvature_df[curvature_df["driver_id"] == driver_id]
        curvature_mean = np.mean(curvature_from_driver_id["curvature"])
        curvature_std = np.std(curvature_from_driver_id["curvature"])
        curvature_summary_rows.append({
            "driver_id": driver_id,
            "curvature_mean": curvature_mean,
            "curvature_std": curvature_std
        })
    return pd.DataFrame(curvature_summary_rows)
