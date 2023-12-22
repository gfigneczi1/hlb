"""Python module to run the whole analysis process for tuning the weights of the speedup formula"""
import logging
import os
import json
import pandas as pd
import numpy as np
import seaborn as sns
from typing import List
from collections import defaultdict
import matplotlib.pyplot as plt
from acceleration_learning.analysis.weighted_formula.weights_solver import WeightsSolver
from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.preprocessing.situation_extraction.time_window_extractor import TimeWindowExtractor
from acceleration_learning.configuration import Configuration
from acceleration_learning.path_configuration import PathConfiguration, TripFile
from acceleration_learning.subject_identifier import SubjectIdentifier

def load_speedup_windows(path_configuration: PathConfiguration, n_trips: int) -> List[driving_situation_window]:
    """Load the stored speedup windows and return then"""
    trip_to_speedups = defaultdict(list)
    for extracted_windows_per_trip_folder in os.listdir(path_configuration.preprocessing_extraction_output_directory_path()):
        if len(trip_to_speedups) == n_trips - 1:
            break
        trip_file = TripFile(extracted_windows_per_trip_folder)
        extractor = TimeWindowExtractor(trip_file.file_name_without_extension)
        windows = extractor.get_saved_windows(
            windows_information_file_path=path_configuration.preprocessing_extraction_information_output_file_path(trip_file),
            input_time_window_directory_path=path_configuration.trip_file_extraction_output_directory_path(trip_file)
        )
        if len(windows[DrivingEventEnum.SPEED_UP]) == 0:
            continue
        trip_to_speedups[trip_file.file_name_without_extension] = windows[DrivingEventEnum.SPEED_UP]
    return trip_to_speedups

def prepare_formula_parameters(
    trip_to_speedups: dict,
    weights_solver: WeightsSolver,
    signal_profile: dict,
    subject_identifier: SubjectIdentifier,
    output_parameters_by_driver_file_path: str
    ) -> None:
    """Compute and save the formula parameters from the given speedup_windows"""
    trip_to_list_of_variables = defaultdict(list)
    for trip_name, speedups in trip_to_speedups.items():
        trip_to_list_of_variables[trip_name] = weights_solver.generate_features_dataset(
            speedup_windows=speedups,
            signal_profile=signal_profile
        )
    subject_to_variables = defaultdict(list)
    for exported_avt_data_trip_name_folder, variables in trip_to_list_of_variables.items():
        trip_name = '_'.join(exported_avt_data_trip_name_folder.split("_")[:3])
        logging.debug("Converting trip: %s", trip_name)
        subject = subject_identifier.get_driver_id(trip_name)
        subject_to_variables[subject] += variables
    logging.info("Saving the parameters with subject into %s", output_parameters_by_driver_file_path)
    with open(output_parameters_by_driver_file_path, "w", encoding='utf-8') as subject_variables_list_file:
        json.dump(subject_to_variables, subject_variables_list_file)

def weighted_formula_tuning_logic(configuration: Configuration) -> None:
    """Tune the weights of the weighted formula of the aggressiveness.
    Generate intermediate file and result file in the dedicated analysis folder.
    Assume that the speedup windows are already stored in the preprocessing output folder
    """
    logging.info("Configuring and setting up.")
    path_configuration = configuration.path_configuration
    subject_identifier = SubjectIdentifier(path_configuration.driver_trip_mapping_file_path())
    subject_variables_list_file_path = os.path.join(
        path_configuration.analysis_weighted_formula_output_directory_path(),
        "formula_parameters_per_subject.json"
    )
    weight_output_file_path = os.path.join(
        path_configuration.analysis_weighted_formula_output_directory_path(),
        "formula_weights.json"
    )
    trip_to_speedups = load_speedup_windows(
        path_configuration=path_configuration,
        n_trips=100
    )
    # Tune weights
    logging.info("Computing formula's parameters")
    solver = WeightsSolver()
    prepare_formula_parameters(
        trip_to_speedups=trip_to_speedups,
        weights_solver=solver,
        signal_profile=configuration.signal_profile,
        subject_identifier=subject_identifier,
        output_parameters_by_driver_file_path=subject_variables_list_file_path
    )
    # Optimize
    logging.info("Tuning the formula's weights")
    weights = solver.find_optimum(subject_variables_list_file_path)
    logging.info("Saving weights into %s", weight_output_file_path)
    solver.save_weights(weight_output_file_path, weights, DrivingEventEnum.SPEED_UP)
    # Display driver distribution
    logging.info("Displaying the aggressiveness's distribution among drivers")
    variables_per_driver = dict()
    with open(subject_variables_list_file_path, "r", encoding='utf-8') as json_file:
        variables_per_driver = json.load(json_file)
    loads_per_driver = defaultdict(list)
    for driver_name, trip_variables_list in variables_per_driver.items():
        variables_means = []
        variables_list = np.transpose(np.array(trip_variables_list))
        # variables_list: [
        #   [j1a, j1b, ...],
        #   [p1a, p1b, ...],
        #   ...]
        for variables in variables_list:
            variables_means.append(np.mean(variables))
        # variables_means: [
        #   j1_mean,
        #   p1_mean,
        #   ...]
        loads_per_driver["driver"].append(driver_name)
        loads_per_driver["aggressiveness"].append(variables_means @ weights)
    df = pd.DataFrame.from_dict(loads_per_driver)
    sns.displot(loads_per_driver, x="aggressiveness", kind="hist")
    plt.xlabel("Aggressiveness")
    plt.title("Distribution of the aggressiveness per driver")
    logging.info('Saving the figure of the distribution')
    plt.savefig(os.path.join(
        path_configuration.analysis_weighted_formula_output_directory_path(),
        'aggressiveness_distribution_among_driver.png'
    ))
    plt.show()
