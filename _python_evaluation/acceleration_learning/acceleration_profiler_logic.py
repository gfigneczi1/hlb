import os
import sys
import time
import logging
from dataclasses import dataclass
import pandas as pd

from acceleration_learning.preprocessing.situation_extraction.time_window_extractor import TimeWindowExtractor
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.preprocessing.driving_event_metadata_analysis import driving_event_reason_analysis
from acceleration_learning.preprocessing.map.tile_id_extraction import extract_tile_id
from acceleration_learning.preprocessing.map.extract_map_data import extract_map_data
from acceleration_learning.modeling.matlab_estimation_helper import matlab_estimation_helper
from acceleration_learning.modeling.acceleration_profile_modeler import acceleration_profile_modeler

from acceleration_learning.path_configuration import PathConfiguration
from acceleration_learning.configuration import Configuration
from acceleration_learning.trip_file import TripFile
from acceleration_learning.subject_identifier import SubjectIdentifier
from acceleration_learning.process_step import ProcessStep
from acceleration_learning.modeling.estimation_step import EstimationStep
from utils.path import build_directory_tree

@dataclass
class TripInformation:
    """Information relative to the trip data container"""
    subject_id: int
    trip_file: TripFile

def required_signals_are_available(signal_profile: dict) -> bool:
    """Return True if the signal required for the "preprocess" step are in the data, False otherwise"""
    if signal_profile["Time"] == "unavailable":
        return False
    if signal_profile["Longitudinal_Velocity"] == "unavailable":
        return False
    if signal_profile["Longitudinal_Acceleration"] == "unavailable":
        return False
    if signal_profile["YawRate"] == "unavailable":
        return False
    return True

def build_file_project_hierarchy(path_configuration: PathConfiguration, trip_file: TripFile):
    """Create the necessary directories to store the steps' outputs"""
    directories_paths = list(filter(lambda path: os.path.isdir(path), path_configuration.path_configuration.values()))
    one_trip_specific_paths = [
        os.path.join(
            path_configuration.trip_file_extraction_output_directory_path(trip_file)
        ),
        os.path.join(path_configuration.preprocessing_map_output_directory_path(), trip_file.file_name)
    ]
    paths_to_build = directories_paths + one_trip_specific_paths
    build_directory_tree(paths_to_build)
    logging.info("The following paths are existing :%s", ",".join(paths_to_build))

def extraction(csv_trip_file: TripFile, signal_profile: dict, subject_id: str, configuration: Configuration) -> dict:
    """Perform the situation extraction
    Return a dictionary of the extractions situations per situation
    """
    trip_data_df = pd.read_csv(csv_trip_file.file_path)
    time_window_extractor = TimeWindowExtractor(csv_trip_file.file_name)
    situation_windows = time_window_extractor.makeWindows(trip_data_df, signal_profile)
    overall_window_index = 0
    window_dict_list = []
    for situation, windows in situation_windows.items():
        logging.info("%d %s situations have been extracted.", len(windows), situation.value)
        for i, window in enumerate(windows):
            output_file = configuration.trip_file_extraction_output_directory_path(csv_trip_file)
            time_window_extractor.save_window_data(
                window=window,
                index=i,
                output_window_directory_path=output_file
            )
            window_dict_list.append(time_window_extractor.make_window_dict(
                window,
                i,
                overall_window_index,
                subject_id
            ))
            overall_window_index += 1
    time_window_extractor.save_window_info(
        window_dict_list=window_dict_list,
        window_information_file_path=configuration.preprocessing_extraction_information_output_file_path(csv_trip_file)
    )
    return situation_windows

def map_signals_are_available(signal_profile: dict) -> bool:
    """Return True if the signal required for the "map" step are in the data, False otherwise"""
    if signal_profile["Longitude"] == "unavailable":
        return False
    if signal_profile["Latitude"] == "unavailable":
        return False
    return True

def modeling(
    configuration: Configuration,
    situation_windows: dict,
    input_windows_directory_path: str,
    trip_information: TripInformation,
    time_window_extractor: TimeWindowExtractor,
    ):
    """Perform speed up modeling and parameters estimation"""
    matlab_helper = matlab_estimation_helper()
    estimation_settings = configuration.estimation_settings
    path_configuration = configuration.path_configuration
    if estimation_settings["Estimation_Step"] == EstimationStep.PYTHON_ESTIMATION.value.lower():
        if len(situation_windows[DrivingEventEnum.SPEED_UP]):
            #modeling parameters
            accel_profile_modeler = acceleration_profile_modeler()
            trip_parameters, all_speedup_parameters = accel_profile_modeler.trip_parameter_estimation(
                input_windows_directory_path=input_windows_directory_path,
                signal_profile=configuration.signal_profile,
                model=estimation_settings["model"]
            )
            #updating the windows
            overall_window_index = 0
            speedup_counter = 0
            window_dict_list = []
            if all_speedup_parameters.shape[0] > 0:
                for windows in situation_windows.values():
                    for i, window in enumerate(windows):
                        if window.classification.value == DrivingEventEnum.SPEED_UP.value:
                            window.parameters = list(all_speedup_parameters[speedup_counter, :])
                            speedup_counter += 1
                        window_dict_list.append(time_window_extractor.make_window_dict(
                            window=window,
                            index=i,
                            window_index=overall_window_index,
                            subject_id=trip_information.subject_id
                        ))
                        overall_window_index += 1
                time_window_extractor.save_window_info(
                    window_dict_list=window_dict_list,
                    window_information_file_path=path_configuration.preprocessing_extraction_information_output_file_path(
                        trip_file=trip_information.trip_file
                    )
                )
                accel_profile_modeler.save_trip_parameter_output(
                    output_parameter_file_path=path_configuration.modeling_output_trip_parameter_path(
                        trip_file=trip_information.trip_file
                    ),
                    trip_parameters=trip_parameters,
                    model=estimation_settings["model"]
                )
                logging.info(
                    "The extracted average parameters for speedups have been saved in the '%s' file",
                    path_configuration.modeling_output_trip_parameter_path(
                        trip_file=trip_information.trip_file
                    )
                )
            else:
                logging.info("No parameters could be estimated.")
    elif estimation_settings["Estimation_Step"] == EstimationStep.MATLAB_ESTIMATION_PREPARATION.value.lower():
        matlab_helper.prepare_matlab_estimation(
            mapping_file_path=configuration.subject_trip_mapping_file_path,
            input_multiple_trip_windows_directory_path=path_configuration.preprocessing_extraction_output_directory_path(),
            input_output_postfixed_trip_information_file_path=path_configuration.modeling_trip_information_input_file_path(
                postfix_name=estimation_settings["postfix"]
            )
        )
    elif estimation_settings["Estimation_Step"] == EstimationStep.MATLAB_COMMANDS_PRINT.value.lower():
        matlab_helper.print_matlab_commands(
            estimation_settings=estimation_settings,
            signal_profile=configuration.signal_profile,
            postfixed_matlab_model_parameters_file_path=path_configuration.modeling_matlab_model_parameters_input_file_path(
                configuration.estimation_settings["postfix"]
            ),
            input_postfixed_trip_information_file_path=path_configuration.modeling_trip_information_input_file_path(
                configuration.estimation_settings["postfix"]
            ),
            extraction_output_directory_path=path_configuration.preprocessing_extraction_output_directory_path()
        )
    elif estimation_settings["Estimation_Step"] == EstimationStep.MATLAB_RESULTS_SAVE.value.lower():
        matlab_helper.save_matlab_estimates_to_info_files(
            estimation_settings=estimation_settings,
            postfixed_matlab_model_parameters_file_path=path_configuration.modeling_matlab_model_parameters_input_file_path(
                postfix_name=estimation_settings["postfix"]
            ),
            postfixed_trip_information_file_path= path_configuration.modeling_trip_information_input_file_path(
                postfix_name=estimation_settings["postfix"]
            ),
            window_information_file_path=path_configuration.preprocessing_extraction_information_output_file_path(
                trip_file=trip_information.trip_file
            ),
            trip_parameter_output_file_path=path_configuration.modeling_output_trip_parameter_path(
                trip_file=trip_information.trip_file
            )
        )

def process_step_help()-> str:
    """Return an help message explaining the step of the process_trip_file function"""
    help_message = " step options\n"
    help_message += "- temp_folder_preparation: if a directory needed in the running of the remaining script and doesn't exist, create it\n"
    help_message += "It is recommended to delete the extracted data for a trip or to move the data, when a new trip is input, \n"
    help_message += "since otherwise the files from the old trip will be overwritten. Recommendation: move the entire folder and run this step to generate a new one.\n"
    help_message += "- extraction: detection of driving situations and extraction of time sequences (so called windows)\n"
    help_message += "    -> saves csv with data of window sequence\n"
    help_message += "    -> saves json with meta information about window sequences\n"
    help_message += "- tile_id_extraction\n"
    help_message += "    -> prints out IDs of all tiles of interest for the map export in the terminal\n"
    help_message += "- map_extraction\n"
    help_message += "    -> saves csv with the original window sequence data and the additionally matched map information\n"
    help_message += "- base_analysis\n"
    help_message += "    -> adds the reasons to the metadata json about the sequences\n"
    help_message += "- map_analysis\n"
    help_message += "    -> adds the reasons to the metadata json about the sequences including map related reasons\n"
    help_message += "- modeling\n"
    help_message += "    -> adds parameters to the metadata json about sequences for speedups\n"
    help_message += "    -> saves a file with the average speedup parameters of the trip\n"
    help_message += "- preprocessing: _temp_folder_preparation + preprocessing + base_analysis\n"
    return help_message

def process_trip_file(input_file_path: str, step: str, configuration: Configuration):
    """Run the data trip analysis pipeline for one trip file
    It is assumed that the input file exist"""
    signal_profile = configuration.signal_profile
    path_configuration = configuration.path_configuration
    use_map = False
    trip_file = TripFile(file_path=input_file_path)
    subject_id = SubjectIdentifier.get_subject_id(
        trip_name=trip_file.file_name,
        mapping_file_path=configuration.subject_trip_mapping_file_path
    )
    time_window_extractor = TimeWindowExtractor(trip_file.file_name_without_extension)
    reason_analyzer = driving_event_reason_analysis()
    logging.info("Process of %s starts.", trip_file.file_name_without_extension)
    if not required_signals_are_available(signal_profile):
        logging.error("The signal profile does not have the required signals: Time, Longitudinal_Velocity, Longitudinal_Acceleration, YawRate")
        sys.exit()

    #preparing the temp folder in case it is empty
    build_file_project_hierarchy(path_configuration=path_configuration, trip_file=trip_file)
    logging.info("Project directory hierarchy is built.")

    ### Extraction part ###
    #extraction: creating the window sequence elements
    situation_windows =  {}
    if step in  [ProcessStep.TIME_WINDOW_EXTRACTION,  ProcessStep.PREPROCESSING]:
        logging.info("Time windows extractions starts.")
        start_time = time.time()
        situation_windows = extraction(
            csv_trip_file=trip_file,
            signal_profile=signal_profile,
            subject_id=subject_id,
            configuration=path_configuration
        )
        end_time = time.time()
        time_spent = end_time - start_time
        logging.info("Time windows extractions ended (time spent: %d seconds).", time_spent)

    ### MAP USAGE ###
    # It is not recommended to do map stuff with the folder but only with the individual file since
    # some copy commands may get forgotten
    if map_signals_are_available(signal_profile):
        # Tile ID extraction for map export
        if step in [ProcessStep.TITLE_ID_EXTRACTION]:
            logging.info("Map tiles extraction starts.")
            extract_tile_id(
                input_trip_folder_path=path_configuration.trip_file_extraction_output_directory_path(
                    trip_file=trip_file
                ),
                output_geographic_data_path=path_configuration.geographic_data_output_file_path(
                    trip_file=trip_file
                ),
                signal_profile=signal_profile
            )
            logging.info("Map tiles extraction ended.")

        # Extract information to the map, match it to the original data sequence and fuse it to
        # one csv file
        if step in [ProcessStep.MAP_EXTRACTION]:
            logging.info("Map data extraction starts")
            extract_map_data(
                geographic_file_path=path_configuration.geographic_data_output_file_path(),
                map_info_file_path=path_configuration.map_info_path(),
                input_trip_data_folder_path=path_configuration.input_directory_path(),
                output_map_data_folder_path=path_configuration.geographic_data_output_directory_path(),
                signal_profile=signal_profile
            )
            logging.info("Map data extraction ended.")

    ### ANALYSIS PART ###
    input_windows_directory_path = path_configuration.trip_file_extraction_output_directory_path(trip_file=trip_file)
    #analyze the information to specify the reason for the event
    if step in [ProcessStep.MAP_ANALYSIS] and map_signals_are_available(signal_profile):
        use_map = True
        input_windows_directory_path = path_configuration.geographic_trip_directory_output_directory_path(trip_file=trip_file)
    
    # Load the windows back in if only the latter stages are executed
    if step in [ProcessStep.MAP_ANALYSIS, ProcessStep.BASE_ANALYSIS, ProcessStep.MODELING]:
        logging.info("Loading the time windows starts")
        start_time = time.time()
        situation_windows = time_window_extractor.get_saved_windows(
            windows_information_file_path=path_configuration.preprocessing_extraction_information_output_file_path(
                trip_file=trip_file
            ),
            input_time_window_directory_path=input_windows_directory_path
        )
        end_time = time.time()
        time_spent = end_time - start_time
        logging.info("Loading the time windows ended (it took %f seconds).", time_spent)

    #calculate the Level of Acceleration Dynamics value
    if step in [ProcessStep.BASE_ANALYSIS, ProcessStep.MAP_ANALYSIS, ProcessStep.PREPROCESSING]:
        logging.info("Metadata analysis starts.")
        overall_window_index = 0
        window_dict_list = []
        for windows in situation_windows.values():
            for i, window in enumerate(windows):
                #get reasons for driving event
                window = reason_analyzer.execute(window, signal_profile, use_map=use_map)
                #update the info file
                window_dict_list.append(time_window_extractor.make_window_dict(
                    window=window,
                    index=i,
                    window_index=overall_window_index,
                    subject_id=subject_id
                ))
                overall_window_index += 1
        time_window_extractor.save_window_info(
            window_dict_list=window_dict_list,
            window_information_file_path=path_configuration.preprocessing_extraction_information_output_file_path(
                trip_file=trip_file
            )
        )
        logging.info("Metadata analysis ended.")

    #the final modeling step to extract parameters from creating the acceleration targets
    if step in [ProcessStep.MODELING]:
        logging.info("Modeling starts.")
        modeling(
            configuration=configuration,
            situation_windows=situation_windows,
            input_windows_directory_path=input_windows_directory_path,
            trip_information=TripInformation(subject_id=subject_id, trip_file=trip_file),
            time_window_extractor=time_window_extractor
        )
        logging.info("Modeling ended.")
    logging.info("The %s data file has been processed sucessfully", trip_file.file_name)
    logging.info("Congratulations! You have reached the end!")
