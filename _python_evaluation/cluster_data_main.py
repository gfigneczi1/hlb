import os
import pathlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
from collections import defaultdict
import statistics
from sklearn import preprocessing
# Performing clustering
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.cluster import DBSCAN

# Imports different color maps for plots
import matplotlib.cm as cm

from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.analysis.clustering.cluster_methods import cluster_with_kmeans, cluster_with_gaussian_mixture_model, distance_time_wrapping, cluster_with_time_series_kmeans, cluster_with_kmedoid
from utils.statistics import filter_outliers_all_columns

TEMP_FOLDER_PATH = os.path.join(r'C:\git\kdp_hlb_evalframework\_temp')
CLUSTERING_FOLDER_PATH = os.path.join(r'C:\git\kdp_hlb_evalframework\_temp\analysis\clustering')

def min_max_norm(value, min, max):
    norm = (value - min) / (max - min)
    return norm

def plot_2D(df, columns):
    first_prop = columns[0]
    second_prop = columns[1]
    plt.title(f"{first_prop} vs {second_prop}")
    plt.scatter(df[first_prop], df[second_prop])
    plt.xlabel(first_prop)
    plt.ylabel(second_prop)
    plt.show()

def plot_hist(df, column):
    plt.hist(df[column], bins="auto")
    plt.ylabel("Number of samples")
    plt.xlabel(column)
    plt.title("Distribution of samples of " + column)
    plt.show()

def elbow_method(df, max_iteration):
    distortions = []
    K = range(1, max_iteration+1)
    for k in K:
        model = KMeans(n_clusters=k).fit(df)
        distortions.append(model.inertia_)
    plt.plot(K, distortions, 'go-')
    plt.xlabel("Number of clusters")
    plt.ylabel("Distortion")
    plt.title('The Elbow Method for cluster number')
    plt.show()

def plot_2D_kmeans_cluster(df, columns, model, centers):
    print(centers)
    predicted = model.predict(df)
    first_col = columns[0]
    second_col = columns[1]
    plt.title(f"Comparison of the parameters: {first_col} vs {second_col}")
    plt.scatter(df[first_col], df[second_col], c=predicted, cmap=cm.rainbow)
    plt.scatter(centers[:, df.columns.get_loc(first_col)], centers[:, df.columns.get_loc(second_col)], c='black', s=200, alpha=0.5)
    plt.xlabel(first_col)
    plt.ylabel(second_col)
    plt.show()

def normalize_std(df) -> pd.DataFrame:
    df_std = preprocessing.StandardScaler().fit(df).transform(df)
    df_std = pd.DataFrame(df_std)
    df_std.columns = df.columns
    return df_std

def compute_aggressiveness_parameter(speed_up_window: driving_situation_window) -> list:
    data: pd.DataFrame = speed_up_window.data
    if data.shape[0] == 0:
        raise ValueError("dataframe empty")
    # Time to peak
    peak_instant_index = data["LongitudinalVelocity_0x344_20"].idxmax()
    peak_time = data["ts_s_0"].iloc[peak_instant_index] - data["ts_s_0"].iloc[0]
    # Time to settle
    settling_time_in_s = data["ts_s_0"].iloc[-1] - data["ts_s_0"].iloc[peak_instant_index]
    # Velocity overshoot
    initial_velocity = data["LongitudinalVelocity_0x344_20"].iloc[0]
    final_velocity = data["LongitudinalVelocity_0x344_20"].iloc[-1]
    velocity_change = final_velocity - initial_velocity
    maximal_velocity = max(data["LongitudinalVelocity_0x344_20"])
    velocity_overshoot = maximal_velocity - final_velocity
    velocity_overshoot_ratio = velocity_overshoot / velocity_change
    # Level of acceleration
    maximal_acceleration = max(data["ActualAcceleration_0x17D_100"])
    # "Kneel"
    peak_velocity = data["LongitudinalVelocity_0x344_20"].iloc[peak_instant_index]
    start_velocity = data["LongitudinalVelocity_0x344_20"].iloc[0]
    n_samples = peak_instant_index + 1
    straight_line_from_start_to_peak = np.linspace(start=start_velocity, stop=peak_velocity, num=n_samples)
    velocity_from_start_to_peak = data["LongitudinalVelocity_0x344_20"][0:n_samples]
    diff_sum = 0
    for straight_line_velo, real_velo in zip(straight_line_from_start_to_peak, velocity_from_start_to_peak):
        diff_sum += real_velo - straight_line_velo
    kneel_integral = diff_sum / n_samples
    # Jerk
    return settling_time_in_s, velocity_overshoot_ratio, maximal_acceleration, peak_time, kneel_integral

def get_zeta_from_overshoot(os_target):
    n = 1000
    for i in range(0, n):
        zeta = i / n
        os = np.exp(-np.pi*zeta/np.sqrt(1-zeta**2))
        if os_target - 0.001<= os and os <= os_target + 0.001:
            return zeta
    return 0

def extracted_basic_properties_and_save_json():
    windows_data_folder_path = os.path.join(TEMP_FOLDER_PATH, "preprocessing", "extraction", "window_data")
    # Extract windows
    parameters_per_trip = []
    trips = dict()
    for i, trip_directory in enumerate(os.listdir(windows_data_folder_path)):
        print(f"i:{i}")
        # Load saved speedup windows
        speed_ups = []
        for window_file_name in os.listdir(os.path.join(windows_data_folder_path, trip_directory)):
            if "speed_up" in window_file_name:
                data_path = str(os.path.join(windows_data_folder_path, trip_directory, window_file_name))
                data = pd.read_csv(data_path)
                if data.shape[0] == 0:
                    continue
                window = driving_situation_window(0, data.shape[0]-1, data, DrivingEventEnum.SPEED_UP)
                speed_ups.append(window)
        # Compute aggressiveness per trip
        if len(speed_ups) == 0:
            continue
        settle_times = []
        overshoots = []
        accelerations = []
        peaks = []
        kneels = []
        for speed_up in speed_ups:
            settle_time, overshoot, acceleration, peak_time, kneel = compute_aggressiveness_parameter(speed_up)
            settle_times.append(settle_time)
            peaks.append(peak_time)
            overshoots.append(overshoot)
            accelerations.append(acceleration)
            kneels.append(kneel)
        parameters = [np.mean(peak_time), np.mean(overshoots), np.mean(accelerations), np.mean(settle_times), np.mean(kneels)]
        parameters_per_trip.append(parameters)
        trips[trip_directory] = parameters
    # Group by driver
    database_filename = "speedup_parameters_and_properties.csv"
    database_path = str(os.path.join(TEMP_FOLDER_PATH, database_filename))
    database_speedup_parameters_analysis: pd.DataFrame = pd.read_csv(database_path)
    trip_per_driver = dict()
    for trip_name, properties in trips.items():
        current_trip_row = database_speedup_parameters_analysis.loc[database_speedup_parameters_analysis["trip_name"] == trip_name]
        if current_trip_row.empty:
            continue
        current_driver = current_trip_row["driver_id"].iloc[0]
        if not current_driver in trip_per_driver:
            trip_per_driver[current_driver] = []
        trip_per_driver[current_driver].append(properties)
    # Save parameters
    with open(os.path.join(CLUSTERING_FOLDER_PATH, "parameters.json"), "w") as param_file:
        json.dump(parameters_per_trip, param_file)
    # Save parameters per trip
    with open(os.path.join(CLUSTERING_FOLDER_PATH, "parameters_per_trip.json"), "w") as param_file:
        json.dump(trips, param_file)
    # Save
    with open(os.path.join(CLUSTERING_FOLDER_PATH, "parameters_list_per_driver.json"), "w") as param_file:
        json.dump(trip_per_driver, param_file)

def main_for_cluster_on_basic_properties():
    """ Try a light-weighted approach to cluster speed ups on curve properties (acceleration, peak time, settle time) """
    #extracted_basic_properties_and_save_json()
    # Load parameters
    parameters_per_trip = []
    with open(os.path.join(CLUSTERING_FOLDER_PATH, "parameters.json"), "r", encoding='utf-8') as param_file:
        parameters_per_trip = json.load(param_file)
    peak_times = []
    overshoots = []
    accelerations = []
    settle_times = []
    for parameters in parameters_per_trip:
        peak_time, overshoot, acceleration, settle_time, kneel = parameters
        peak_times.append(peak_time)
        overshoots.append(overshoot)
        accelerations.append(acceleration)
        settle_times.append(settle_time)
    
    d = {
        "peak_time" : peak_times,
        "overshoot" : overshoots,
        "acceleration" : accelerations,
        "settle_time" : settle_times
    }
    data = pd.DataFrame(d)
    
    # Remove outliers
    def filter(df, a):
        for col in df:
            quantile_low = df[col].quantile(a)
            quantile_high  = df[col].quantile(1-a)

            df = df[(df[col] < quantile_high) & (df[col] > quantile_low)]
        return df
    
    data = filter(data, 0.01)
        
    peak_times = []
    overshoots = []
    accelerations = []
    settle_times = []
    parameters_per_trip = data.to_numpy()
    for parameters in parameters_per_trip:
        peak_time, overshoot, acceleration, settle_time = parameters
        peak_times.append(peak_time)
        overshoots.append(overshoot)
        accelerations.append(acceleration)
        settle_times.append(settle_time)
        
    # Normalization
    normalized_parameters = []
    min_peak = min(peak_times)
    max_peak = max(peak_times)
    min_acc = min(accelerations)
    max_acc = max(accelerations)
    min_settle = min(settle_times)
    max_settle = max(settle_times)
    for parameters in parameters_per_trip:
        peak_time, overshoot, acceleration, settle_time = parameters
        norm_peak_time = min_max_norm(peak_time, min_peak, max_peak)
        norm_accel = min_max_norm(acceleration, min_acc, max_acc)
        norm_settle_time = min_max_norm(settle_time, min_settle, max_settle)
        normalized_parameters.append([norm_settle_time, overshoot, norm_accel, norm_peak_time])
    
    # Remove time to settle
    for row in normalized_parameters:
        del row[-1]

    # Clustering
    medoids = cluster_with_kmedoid(normalized_parameters, 3)
    
    # Step response visualization  
    # Formula from https://jckantor.github.io/CBE30338/03.06-Second-Order-Models.html
    def underdamped(K, tau, zeta, period):
        t = np.linspace(0,period)
        c = np.cos(np.sqrt(1-zeta**2)*t/tau)
        s = np.sin(np.sqrt(1-zeta**2)*t/tau)
        
        y = K*(1 - np.exp(-zeta*t/tau)*(c + zeta*s/np.sqrt(1-zeta**2)))
        return t, y
    
    plt.figure(1)
    for i, medoid in enumerate(medoids):
        tp = medoid[0]
        overshoot = medoid[1]
        acc = medoid[2]
        
        peak_time = tp
        zeta = get_zeta_from_overshoot(overshoot)
        tau = peak_time / np.pi * np.sqrt(1-zeta**2)

        period = peak_time*2
        x, y = underdamped(1, tau, zeta, period)
        color = ""
        if i == 0:
            color = "red"
        elif i == 1:
            color = "blue"
        else:
            color = "green"
        plt.plot(x, y, color=color)
    plt.grid()
    plt.show()

def main_for_analyze_model_parameters():
    """ Cluster step by step the parameters database of Samuel"""
    # Load the data
    homedir = str(pathlib.Path(__file__).resolve().parent.parent)
    avt_data_folder_path = os.path.join(homedir, "_temp")
    database_filename = "filtered_parameters.csv"
    database_path = str(os.path.join(homedir, "_temp", database_filename))
    database = pd.read_csv(database_path)
    # Focus on one category of speed up (for example, zero to middle speed)
    velocity_classes = [
        "zero_to_low_speed_up",     # 0
        "zero_to_mid_speed_up",     # 1
        "zero_to_high_speed_up",    # 2
        "low_to_low_speed_up",      # 3
        "low_to_mid_speed_up",      # 4
        "low_to_high_speed_up",     # 5
        "mid_to_mid_speed_up",      # 6
        "mid_to_high_speed_up",     # 7
        "high_to_high_speed_up"     # 8
    ]
    def get_cmap(n, name='hsv'):
        '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
        RGB color; the keyword argument name must be a standard mpl colormap name.'''
        return plt.cm.get_cmap(name, n)
    speedup_parameters_to_plot = ["alfa", "OS", "Ts"]
    minimum_number_of_speedups = 30
    n_drivers_target = 5
    for velocity_class in velocity_classes:
        database_where_velocity_class_is_zero_to_mid = database.loc[database.velocity_class == velocity_class]
        print("database_where_velocity_class_is_zero_to_mid", len(database_where_velocity_class_is_zero_to_mid))
        # Sort the trip by driver
        ## Filter out the drivers with little data
        ## Filter to only have n_drivers_target drivers
        drivers_having_minimum_of_recorded_speedup = database_where_velocity_class_is_zero_to_mid[
            database_where_velocity_class_is_zero_to_mid["driver_id"]
                .map(database_where_velocity_class_is_zero_to_mid["driver_id"].value_counts()) > minimum_number_of_speedups
        ]
        while len(set(drivers_having_minimum_of_recorded_speedup["driver_id"])) > n_drivers_target:
            minimum_number_of_speedups += 5
            drivers_having_minimum_of_recorded_speedup = database_where_velocity_class_is_zero_to_mid[
            database_where_velocity_class_is_zero_to_mid["driver_id"]
                .map(database_where_velocity_class_is_zero_to_mid["driver_id"].value_counts()) > minimum_number_of_speedups
            ]
        if len(set(drivers_having_minimum_of_recorded_speedup["driver_id"])) < n_drivers_target:
            minimum_number_of_speedups -= 5
            drivers_having_minimum_of_recorded_speedup = database_where_velocity_class_is_zero_to_mid[
            database_where_velocity_class_is_zero_to_mid["driver_id"]
                .map(database_where_velocity_class_is_zero_to_mid["driver_id"].value_counts()) > minimum_number_of_speedups
            ]
        print("minimum_number_of_speedups found", minimum_number_of_speedups)
        print("drivers_having_minimum_of_recorded_speedup", len(drivers_having_minimum_of_recorded_speedup["driver_id"]))
        # Plot the scatter of the drivers' speedup
        ax = plt.axes(projection ="3d")
        first_param = speedup_parameters_to_plot[0]
        second_param = speedup_parameters_to_plot[1]
        third_param = speedup_parameters_to_plot[2]
        drivers = set(drivers_having_minimum_of_recorded_speedup["driver_id"])
        n_drivers = len(drivers)
        print("n_driver", n_drivers)
        cmap = get_cmap(n_drivers)
        for driver_i, driver in enumerate(drivers):
            x = drivers_having_minimum_of_recorded_speedup.loc[drivers_having_minimum_of_recorded_speedup.driver_id == driver][first_param]
            y = drivers_having_minimum_of_recorded_speedup.loc[drivers_having_minimum_of_recorded_speedup.driver_id == driver][second_param]
            z = drivers_having_minimum_of_recorded_speedup.loc[drivers_having_minimum_of_recorded_speedup.driver_id == driver][third_param]
            ax.scatter(x, y, z, marker=".", color=cmap(driver_i), label=driver)
        ax.set_xlabel(first_param)
        ax.set_ylabel(second_param)
        ax.set_zlabel(third_param)
        plt.grid()
        plt.title(f"Scatter of {velocity_class} class, with params: {first_param}, {second_param}, {third_param}, having {n_drivers} driver")
        plt.show()

def main_for_group_driver_of_basic_param_extracted():
    # # Load parameters
    # parameters_per_trip = dict()
    # with open("parameters_per_trip.json", "r") as param_file:
    #     parameters_per_trip = json.load(param_file)
    
    # # Load the database
    # homedir = str(pathlib.Path(__file__).resolve().parent.parent)
    # avt_data_folder_path = os.path.join(homedir, "_temp")
    # database_filename = "speedup_parameters_and_properties.csv"
    # database_path = str(os.path.join(homedir, "_temp", database_filename))
    # database_speedup_parameters_analysis: pd.DataFrame = pd.read_csv(database_path)
    # # Group by driver
    # trip_per_driver = dict()
    # for trip_name, properties in parameters_per_trip.items():
    #     current_trip_row = database_speedup_parameters_analysis.loc[database_speedup_parameters_analysis["trip_name"] == trip_name]
    #     if current_trip_row.empty:
    #         continue
    #     current_driver = current_trip_row["driver_id"].iloc[0]
    #     if not current_driver in trip_per_driver:
    #         trip_per_driver[current_driver] = []
    #     trip_per_driver[current_driver].append(properties)
    
    # # Save
    # with open("parameters_list_per_driver.json", "w") as param_file:
    #     json.dump(trip_per_driver, param_file)
    
    trip_per_driver = dict()
    with open(os.path.join(CLUSTERING_FOLDER_PATH, "parameters_list_per_driver.json"), "r") as param_file:
        trip_per_driver = json.load(param_file)
    fig = plt.figure(0)
    ax = plt.axes(projection ="3d")
    # Plot the scatter of the drivers' speedup
    def get_cmap(n, name='hsv'):
        '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
        RGB color; the keyword argument name must be a standard mpl colormap name.'''
        return plt.cm.get_cmap(name, n)
    cmap = get_cmap(len(trip_per_driver))
    averages = defaultdict(list)
    drivers = []
    for driver_i, (driver, parameters_list) in enumerate(trip_per_driver.items()):
        drivers.append(driver)
        peak_times = []
        overshoots = []
        accelerations = []
        settle_times = []
        kneels = []
        for parameters in parameters_list:
            peak_time, overshoot, acceleration, settle_time, kneel = parameters
            peak_times.append(peak_time)
            overshoots.append(overshoot)
            accelerations.append(acceleration)
            settle_times.append(settle_time)
            kneels.append(kneel)
        avg_pt = np.mean(peak_times)
        avg_os = np.mean(overshoots)
        avg_st = np.mean(settle_times)
        avg_ac = np.mean(accelerations)
        avg_kn = np.mean(kneels)
        averages["peak_time"].append(avg_pt)
        averages["overshoot"].append(avg_os)
        averages["settle_time"].append(avg_st)
        averages["acceleration"].append(avg_ac)
        averages["kneel"].append(avg_kn)
        x = accelerations
        y = overshoots
        z = kneels
        ax.scatter(x, y, z, marker=".", color=cmap(driver_i), s = 100, alpha=0.5, label=driver)
    ax.set_xlabel("Acceleration")
    ax.set_ylabel("Overshoot")
    ax.set_zlabel("Kneel")
    plt.grid()
    plt.title("Scatter of speedup's features colored by driver")
    plt.show()
    
    # Filter outliers
    df = pd.DataFrame(averages)
    df = filter_outliers_all_columns(df, 0.2)
    averages = df.to_numpy()
    # Normalize
    df_norm = pd.DataFrame(averages)
    df_norm.columns = df.columns
    print(df_norm)
    df_std = normalize_std(df_norm)
    print(df_norm)
    # Look at the data distribution
    # for param in df_std.columns:
    #     plot_hist(df_std, param)
    # Look at the data in 2D
    # plot_2D(df_std, ["peak_time","overshoot"])
    # plot_2D(df_std, ["peak_time","acceleration"])
    # plot_2D(df_std, ["kneel","acceleration"])
    # plot_2D(df_std, ["kneel","peak_time"])
    # Elbow method to get the number of clusters
    elbow_method(df_std, 10)
    # 2 clusters
    # Perform K-means
    n_clusters = 2
    #model = KMeans(n_clusters=n_clusters).fit(df_std)
    model = GaussianMixture(n_components=n_clusters, n_init=10).fit(df_std)
    #model = DBSCAN(eps=0.3).fit(df_std)
    #centers = model.cluster_centers_
    centers = model.means_
    plot_2D_kmeans_cluster(df_std, ["peak_time", "overshoot"], model, centers)
    plot_2D_kmeans_cluster(df_std, ["settle_time", "overshoot"], model, centers)
    plot_2D_kmeans_cluster(df_std, ["acceleration", "overshoot"], model, centers)
    plot_2D_kmeans_cluster(df_std, ["kneel", "overshoot"], model, centers)
    plot_2D_kmeans_cluster(df_std, ["kneel", "acceleration"], model, centers)

if __name__ == "__main__":
    # extracted_basic_properties_and_save_json()
    main_for_group_driver_of_basic_param_extracted()
    #main_for_cluster_on_basic_properties()
    # Visualize model parameter
    # main_for_analyze_model_parameters()
