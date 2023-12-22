"""Python main script to analyse the speed up velocity curvature."""
import os
import logging
import pandas as pd
from acceleration_learning.analysis.speedup_curvature.speed_up_curvature_computer import compute_speedup_curvature, summary
from acceleration_learning.analysis.speedup_curvature.speedup_curvature_analysis import speedup_curvature_analysis_prototype_main
from acceleration_learning.configuration import Configuration, default_configuration
from acceleration_learning.path_configuration import PathConfiguration
from utils.math import min_max_normalize_data_frame
from utils.path import build_directory_tree

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    configuration: Configuration = default_configuration()
    path_configuration: PathConfiguration = configuration.path_configuration
    build_directory_tree([
        path_configuration.analysis_curvature_output_directory_path()
    ])
    # Compute curvatures
    curvatures_df: pd.DataFrame = compute_speedup_curvature(
        signal_profile=configuration.signal_profile,
        path_configuration=path_configuration,
    )
    curvature_file_path = os.path.join(
        path_configuration.analysis_curvature_output_directory_path(),
        "driver_speedup_curvature_database.csv"
    )
    curvatures_df.to_csv(curvature_file_path)
    logging.info("Saved curvatures at %s", curvature_file_path)
    # Describe the curvatures
    curvature_summary: pd.DataFrame = summary(curvatures_df)
    curvature_summary_file_path = os.path.join(
        path_configuration.analysis_curvature_output_directory_path(),
        "curvature_summary.csv"
    )
    curvature_summary.to_csv(curvature_summary_file_path)
    logging.info("Saved curvatures summary at %s", curvature_summary_file_path)
    #curvature_summary = pd.read_csv(curvature_summary_file_path)
    # Normalize the curvatures
    scaled_curvature_summary_df = min_max_normalize_data_frame(curvature_summary, "curvature_mean")
    scaled_curvature_summary_file_path = os.path.join(
        path_configuration.analysis_curvature_output_directory_path(),
        "scaled_curvature_summary.csv"
    )
    scaled_curvature_summary_df.to_csv(
        scaled_curvature_summary_file_path,
        columns=curvature_summary.columns
    )
    logging.info("Saved scaled curvatures_df summary at %s", scaled_curvature_summary_file_path)

    speedup_curvature_analysis_prototype_main(curvature_file_path)
