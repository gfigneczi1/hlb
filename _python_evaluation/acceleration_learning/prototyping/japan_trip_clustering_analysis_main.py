import os
import json
import pandas as pd
from acceleration_learning.analysis.clustering.cluster_methods import kmean_cluster, gaussian_mixture_cluster
from utils.path import modeling_directory_path, preprocessing_directory_path, temp_path

if __name__ == "__main__":
    trip_driver_mapping_df = pd.read_csv(
        os.path.join(temp_path(), "JP_extreme_experiment_ID_mapping.csv"),
        sep=";"
    )
    parameters_directory_path = os.path.join(
        modeling_directory_path(),
        "trip_parameters"
    )
    trip_directory_path = os.path.join(preprocessing_directory_path(), "extraction", "window_data")
    database_rows = list()
    for parameter_json_file in os.listdir(parameters_directory_path):
        trip_label = "_".join(parameter_json_file.split("_")[:-1]) # remove the "_parameters.json" part
        subject_id = trip_driver_mapping_df[trip_driver_mapping_df["trip_label"] == trip_label]["subject_id"].iloc[0]
        with open(os.path.join(parameters_directory_path, parameter_json_file), "r") as parameters_file:
            parameters_dict = json.load(parameters_file)
            parameters_dict["trip_label"] = trip_label
            parameters_dict["subject_id"] = subject_id
            database_rows.append(parameters_dict)
    driver_trip_parameters_df = pd.DataFrame(database_rows)
    gaussian_mixture_cluster(
        database_df=driver_trip_parameters_df,
        n_clusters=2,
        database_parameters=["D", "alfa"]
    )
