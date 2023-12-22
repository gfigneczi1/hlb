import os, sys
from typing import Callable
import logging
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../..')))
from utils.path import root_path

def from_tripname_get_driver(database_speedup_parameters_analysis, trip_name):
    """Return the driver id of the trip"""
    row = database_speedup_parameters_analysis.loc[
        database_speedup_parameters_analysis["trip_name"] == trip_name
    ]
    if row.empty:
        return None
    return row["driver_id"].iloc[0]

def speedups_to_driver_to_property(property_function: Callable[[pd.DataFrame], float]) -> dict:
    """Return a dict matching the driver id to their trip properties"""
    # Extract windows
    windows_data_folder_path = os.path.join(root_path(), "_temp", "driving_event_window_data")
    database_path = str(os.path.join(root_path(), "_temp", "speedup_parameters_and_properties.csv"))
    database_speedup_parameters_analysis: pd.DataFrame = pd.read_csv(database_path)
    driver_to_property = dict()    # { "driver_id" : [property1, property2, ..]}
    LIMIT_TIME_IN_S = 20
    FRAME_RATE = 30
    n_directories = len(os.listdir(windows_data_folder_path))
    for i, trip_directory in enumerate(os.listdir(windows_data_folder_path)):
        logging.info(f"Processing % (%/%)", trip_directory, i, n_directories)
        # Load saved speedup windows
        property_list = []
        for window_file_name in os.listdir(os.path.join(windows_data_folder_path, trip_directory)):
            if "speed_up" in window_file_name:
                data_path = str(os.path.join(windows_data_folder_path, trip_directory, window_file_name))
                speed_up_df = pd.read_csv(data_path)
                speedup_time = speed_up_df["ts_s_0"].iloc[-1] - speed_up_df["ts_s_0"].iloc[0]
                if speedup_time > LIMIT_TIME_IN_S * FRAME_RATE:
                    continue
                speedup_property = property_function(speed_up_df)
                property_list.append(speedup_property)
        #filtered_properties = filter_long_speedup(property_list)
        driver = from_tripname_get_driver(database_speedup_parameters_analysis, trip_directory)
        driver_to_property[driver] = property_list
    return driver_to_property
