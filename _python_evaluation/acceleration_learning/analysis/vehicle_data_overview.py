import os, sys
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from utils.inputoutput import save_to_json, load_json
from utils.statistics import plot_hist
from utils.path import root_path
from utils.process import from_tripname_get_driver

def compute_length(data: pd.DataFrame) -> float:
    """Compute the length of the trip"""
    if data.empty:
        raise ValueError("data frame empty")
    return data["ts_s_0"].iloc[-1] - data["ts_s_0"].iloc[0]

def speedups_transform_to_driver_to_speedups() -> dict:
    """Return a dict matching the driver id to their trip properties"""
    # Extract windows
    windows_data_folder_path = os.path.join(root_path(), "_temp", "driving_event_window_data")
    database_path = str(os.path.join(root_path(), "_temp", "speedup_parameters_and_properties.csv"))
    database_speedup_parameters_analysis: pd.DataFrame = pd.read_csv(database_path)
    driver_to_speedups = dict()    # { "driver_id" : [property1, property2, ..]}
    n_directories = len(os.listdir(windows_data_folder_path))
    for i, trip_directory in enumerate(os.listdir(windows_data_folder_path)):
        logging.info(f"Processing % (%/%)", trip_directory, i, n_directories)
        # Load saved speedup windows
        speedups = dict()
        for window_file_name in os.listdir(os.path.join(windows_data_folder_path, trip_directory)):
            if "speed_up" in window_file_name:
                data_path = str(os.path.join(windows_data_folder_path, trip_directory, window_file_name))
                speed_up_df = pd.read_csv(data_path)
                speedups[trip_directory + "-" + window_file_name] = speed_up_df.to_dict()
        driver = from_tripname_get_driver(database_speedup_parameters_analysis, trip_directory)
        driver_to_speedups[driver] = speedups
    return driver_to_speedups

def vehicle_data_overview():
    driver_to_speedups: dict =  speedups_transform_to_driver_to_speedups()
    save_to_json(driver_to_speedups, "driver_to_speedups")
    driver_to_length = load_json("driver_to_length")
    print(driver_to_length)
    # Distribution of all drivers
    average_length = list()
    driver_to_n_length = dict()
    all_length = []
    
    for driver, property_list in driver_to_length.items():
        #plot_hist(property_list, f"speedup length of {driver}")
        mean = np.mean(property_list)
        std = np.std(property_list)
        average_length.append(mean)
        driver_to_n_length[driver] = {
            "n": len(property_list),
            "mean": mean,
            "std": std
        }
        for length in property_list:
            all_length.append(length)
    print(driver_to_n_length)
    print("mean", np.mean(all_length), "std", np.std(all_length))
    plot_hist(all_length, "All lengths")
    plot_hist(average_length, "Average length per driver")
    sns.displot(average_length, bins=10)
    plt.show()
    # All distributions per driver
    rows = list()
    for driver, length_list in driver_to_length.items():
        if driver not in ("null", "nan", "None"):
            for length in length_list:
                rows.append([driver, length])
    df = pd.DataFrame(rows)
    df.columns = ["driver", "length"]
    sns.displot(df, x="length", hue="driver", kind="kde", fill=True, common_norm=False)
    #sns.displot(df, x="length", hue="driver", kind="kde", fill=True, common_norm=False) # normalized density on other variable https://seaborn.pydata.org/tutorial/distributions.html#id1
    plt.show()
