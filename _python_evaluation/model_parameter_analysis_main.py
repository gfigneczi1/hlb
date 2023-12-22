"""Main script to only run global model parameter classification"""
import os
from itertools import combinations
import json
import numpy as np
import pandas as pd
from utils.path import input_directory_path, analysis_directory_path, build_directory_tree
from utils.math import min_max_normalize_data_frame
from acceleration_learning.analysis.classification.profiling import create_profile, summarize_profile
from acceleration_learning.analysis.classification.discriminator import Discriminator
from acceleration_learning.modeling.model_parameter import CAPT2ModelParameter
from acceleration_learning.analysis.classification.velocity_class import VelocityClass

if __name__ == "__main__":    
    VELOCITY_CLASSES = list(VelocityClass)
    MODEL_PARAMETERS = list(CAPT2ModelParameter)
    MODEL_PARAMETERS_OUT_DIRECTORY_NAME = "model_parameters"
    MODEL_PARAMETERS_OUTPUT_DIRECTORY_PATH = os.path.join(
        analysis_directory_path(),
        MODEL_PARAMETERS_OUT_DIRECTORY_NAME
    )
    build_directory_tree([
        MODEL_PARAMETERS_OUTPUT_DIRECTORY_PATH
    ])
    model_parameters_path = os.path.join(input_directory_path(), "filtered_parameters.csv")
    model_parameters_df = pd.read_csv(model_parameters_path)
    drivers = model_parameters_df["driver_id"].unique()
    discriminator = Discriminator(database=model_parameters_df)
    # Profile for every driver
    profiles = []
    profiles_df = []
    profiles_np = []
    for driver in drivers:
        profile = create_profile(
            driver_id=driver,
            database=model_parameters_df,
            discriminator=discriminator
        )
        profiles.append(profile)
        profile_df = pd.DataFrame(profile, columns=MODEL_PARAMETERS)
        profiles_df.append(profile_df)
        profile_np = np.matrix(profile)
        profiles_np.append(profile_np)
    # Summarize the profiles
    profile_summaries = []
    for i, profile in enumerate(profiles):
        profile_summaries.append([drivers[i]] + summarize_profile(profile))
    profile_summaries_df = pd.DataFrame(profile_summaries, columns=["driver_id"] + VELOCITY_CLASSES)
    profile_summaries_df_path: str = os.path.join(
        MODEL_PARAMETERS_OUTPUT_DIRECTORY_PATH,
        "profile_summaries_df.csv"
    )
    profile_summaries_df.to_csv(profile_summaries_df_path)
    # Summarize the classes
    profile_of_classes_summaries = []
    for i, profiles in enumerate(profile_summaries):
        profile_of_classes_summaries.append([drivers[i]] + [np.mean(profiles[1:])])
    global_profile_summaries_df = pd.DataFrame(
        profile_of_classes_summaries,
        columns=["driver_id", "global_summary"]
    )
    global_profile_summaries_df_path: str = os.path.join(
        MODEL_PARAMETERS_OUTPUT_DIRECTORY_PATH,
        "global_profile_summaries_df.csv"
    )
    global_profile_summaries_df.to_csv(global_profile_summaries_df_path)
    # Summarize the classes and the parameters
    profile_of_classes_summaries_and_global = []
    for i, profiles in enumerate(profile_summaries):
        profile_of_classes_summaries_and_global.append(
            [drivers[i]] + profiles[1:] + [np.mean(profiles[1:])]
        )
    driver_id_to_classes_summary_and_global_aggressiveness_df = pd.DataFrame(
        profile_of_classes_summaries_and_global,
        columns=["driver_id"] + VELOCITY_CLASSES + ["global_summary"]
    )
    driver_id_to_classes_summary_and_global_aggressiveness_path: str = os.path.join(
        MODEL_PARAMETERS_OUTPUT_DIRECTORY_PATH,
        "driver_id_to_classes_summary_and_global_aggressiveness.csv"
    )
    driver_id_to_classes_summary_and_global_aggressiveness_df.sort_values(
        "global_summary",
        ascending=True
    ).to_csv(driver_id_to_classes_summary_and_global_aggressiveness_path)
    # Normalize between 0 and 1
    scaled_parameters_summary_df = min_max_normalize_data_frame(global_profile_summaries_df, "global_summary")
    scaled_parameters_summary_file_path = os.path.join(
        analysis_directory_path(),
        "model_parameters",
        "scaled_parameters_summary.csv"
    )
    scaled_parameters_summary_df.to_csv(scaled_parameters_summary_file_path)
