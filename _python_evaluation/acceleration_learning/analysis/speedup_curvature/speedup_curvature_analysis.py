"""Python module to analyze the speedup curvature metric"""
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from utils.statistics import plot_hist, filter_outliers_all_columns

def is_curvature_flat(curvature, tolerance, offset):
    """Return true if the curvature is in the flatness tolerance"""
    return abs(curvature) < tolerance + offset

def speedup_curvature_analysis_prototype_main(curvature_database_file_path):
    """Display statistics about speedup's curvature"""
    curvature_database = pd.read_csv(curvature_database_file_path)
    driver_to_curvature = {} # { "driver_id1" : [curvature1, curvature2, ..]}
    for driver in curvature_database['driver_id']:
        driver_to_curvature[driver] = curvature_database[curvature_database['driver_id'] == driver]['curvature'].to_list()
    all_curvatures = []
    for curvature_list in driver_to_curvature.values():
        all_curvatures += curvature_list
    # Filter outliers
    df = pd.DataFrame(all_curvatures)
    df = filter_outliers_all_columns(df, 0.1)
    all_curvatures = df.to_numpy()
    mean_curv = np.mean(all_curvatures)
    std_curv = np.std(all_curvatures)
    median_curv = np.mean(all_curvatures)
    variance_curv = np.var(all_curvatures)
    logging.info(f"mean_curv {mean_curv}, median_curv {median_curv}, std_curv {std_curv}, variance_curv {variance_curv}")
    # Ratio of up, flat and down per driver
    driver_to_curvature_stats = dict()
    for driver, curvature_list in driver_to_curvature.items():
        curvature_stats = {
            "hump": 0,
            "flat" :0,
            "hollow": 0
        }
        for curvature in curvature_list:
            if is_curvature_flat(curvature, std_curv, mean_curv):
                curvature_stats["flat"] += 1
            else:
                if curvature > 0:
                    curvature_stats["hump"] += 1
                else:
                    curvature_stats["hollow"] += 1
        driver_to_curvature_stats[driver] = curvature_stats
    logging.info(f"driver_to_curvature_stats {driver_to_curvature_stats}")
    driver_curvature_classes_stats = {
        "hump": 0,
        "flat": 0,
        "hollow": 0
    }
    for driver, curvature_stats in driver_to_curvature_stats.items():
        driver_curvature_classes_stats[max(curvature_stats, key=curvature_stats.get)] += 1
    logging.info(f"driver_curvature_classes {driver_curvature_classes_stats}")
    # For all ground_truth_driver_aggressiveness
    avg_curvature = list()
    for driver, curvature_list in driver_to_curvature.items():
        mean = np.mean(curvature_list)
        avg_curvature.append(mean)
    plot_hist(avg_curvature, "average curvature per driver")
    data = list()
    ground_truth_driver_aggressiveness = {
        "aggressive": [
            "2016k_103",
            "2016k_198",
            "2016k_097"
        ],
        "calm": [
            "2016k_131",
            "2016k_209",
            "2016k_161"
        ]
    }
    for driver, curvature_list in driver_to_curvature.items():
        if driver != "null" and driver != "nan":
            for curvature in curvature_list:
                behavior = "unspecified"
                color = "grey"
                if driver in ground_truth_driver_aggressiveness["aggressive"]:
                    behavior = "aggressive"
                    color = "red"
                elif driver in ground_truth_driver_aggressiveness["calm"]:
                    behavior = "calm"
                    color = "blue"
            curvature_mean = np.mean(curvature_list)
            data.append([driver, curvature_mean, behavior, color])
    df = pd.DataFrame(data)
    df.columns = ["driver", "curvature", "behavior", "color"]
    plt.legend(loc='upper right')
    sns.catplot(data=df, x="behavior", y="curvature", kind="box")
    sns.displot(df, x="curvature", kind="kde", hue="behavior", fill=True, common_norm=True)
    plt.show()
