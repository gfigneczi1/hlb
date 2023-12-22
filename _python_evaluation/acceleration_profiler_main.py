"""Main main file"""
import os
import sys
import logging
import json
from datetime import datetime
import argparse
from argparse import RawTextHelpFormatter
import toml
from gooey import GooeyParser, Gooey

from utils.path import signal_profiles_path, build_directory_tree
from acceleration_learning.acceleration_profiler_logic import process_trip_file, process_step_help
from acceleration_learning.configuration import Configuration
from acceleration_learning.trip_file import TripFile
from acceleration_learning.path_configuration import PathConfiguration
from acceleration_learning.process_step import ProcessStep
from acceleration_learning.modeling.estimation_step import EstimationStep
from acceleration_learning.modeling.model import Model

def define_argument_with_cli(default_configuration_file_path: str, signal_profiles: list):
    """Setup the arguments for a command line usage"""
    parser = argparse.ArgumentParser(
        description="Analysis the vehicle data",
        formatter_class=RawTextHelpFormatter
    )
    required_argument_parser = parser.add_argument_group('required arguments')
    required_argument_parser.add_argument(
        "-config",
        "--toml_configuration_file_path",
        dest="toml_configuration_file_path",
        type=str,
        default=default_configuration_file_path,
        help="The path to the TOML configuration file describing the paths to temporary folders"
    )
    process_steps = [step.name.lower() for step in list(ProcessStep)]
    required_argument_parser.add_argument(
        '-s',
        '--step',
        required=True,
        type=str,
        dest="process_step",
        default=ProcessStep.PREPROCESSING.name.lower(),
        choices=process_steps,
        help=process_step_help()
    )
    required_argument_parser.add_argument(
        '-im',
        '--InputMode',
        type=str,
        dest="input_mode",
        default='File',
        choices=['File','Folder'],
        help="The File for one single data file or Folder for multiple data files in one same folder as input."
    )
    required_argument_parser.add_argument(
        '-sp',
        '--SignalProfile',
        dest="signal_profile",
        default="Cadillac_CT6",
        choices=signal_profiles,
        help="The signal labels set that are used in the inputted data file."
    )
    required_argument_parser.add_argument(
        '-mp',
        '--MappingFile',
        required=True,
        dest="mapping_file",
        help='CSV file that maps drivers to trips. It must contain "trip_label" and "subject_id" fields. The separator must be a semi coma ";".'
    )
    single_group = required_argument_parser.add_argument_group(
        "Single Input Options",
        "Choose if only one file is analyzed"
    )
    single_group.add_argument(
        '-if',
        '--InputFile',
        type=str,
        dest="input_file",
        help="The path to the data file to process."
    )
    multi_group = required_argument_parser.add_argument_group(
        "Multi Input Options",
        "Choose if an entire folder is analyzed"
    )
    multi_group.add_argument(
        '-iff',
        '--InputFolderPath',
        type=str,
        dest="input_folder_path",
        help="The path to the folder containing one or severals data file to process."
    )
    multi_group.add_argument(
        '-u',
        '--UseAll',
        action='store_true',
        default=False,
        dest="use_all_files",
        help="If True, it applies step on all files in the folder (default: False)"
    )
    multi_group.add_argument(
        '-fis',
        '--FolderInputStart',
        type=int,
        default=0,
        dest="starting_folder_input",
        help="The start index to select the files in a same folder to process (index starts at zero)."
    )
    multi_group.add_argument(
        '-fie',
        '--FolderInputEnd',
        type=int,
        default=1,
        dest="ending_folder_input",
        help="The end index to select the files in a same folder to process."
    )
    parameter_estimation_group = parser.add_argument_group(
        "Parameter Estimation Options",
        "Choose your options for the parameter estimation. These are only relevant if 'modeling' step is chosen."
    )
    estimations_steps = [estimation.value.lower() for estimation in list(EstimationStep)]
    parameter_estimation_group.add_argument(
        '-e',
        '--Estimation_Step',
        type=str,
        dest="estimation_step",
        choices=estimations_steps,
        help="The step to execute in the model parameter estimation step of the process."
    )
    models = [model.value for model in list(Model)]
    parameter_estimation_group.add_argument(
        '-m',
        '--model',
        type=str,
        default=Model.CONSTANT_ACCELERATION_AND_SECOND_ORDER_STEP_RESPONSE.value,
        dest="model",
        choices=models,
        help="Physical model used to represent acceleration patterns."
    )
    matlab_options = parameter_estimation_group.add_argument_group(
        "MATLAB Estimation options",
        "For 'Print_MATLAB_commands'. Requires 'Print_MATLAB_commands'. MATLAB commands automatically saved in clipboard."
    )
    matlab_options.add_argument(
        '-mc',
        '--MATLAB_command',
        type=str,
        dest="matlab_command",
        choices=["parameter_estimation", "save_parameters"],
        help="The MATLAB step needs to be specified"
    )
    matlab_options.add_argument(
        '-ssi',
        '--Starting_Speedup_Index',
        type=int,
        dest="starting_speedup_index",
        help='Based on "trip_infos.csv" from "Prepare_MATLAB_Estimation"'
    )
    matlab_options.add_argument(
        '-esi',
        '--Ending_Speedup_Index',
        type=int,
        dest="ending_speedup_index",
        help='Based on "trip_infos.csv" from "Prepare_MATLAB_Estimation"'
    )
    matlab_options.add_argument(
        '-spf',
        '--Save_postfix',
        dest="save_postfix",
        help='Postfix with which the file gets saved in the save_parameters command.'
    )
    log_levels = ["DEBUG", "INFO", "WARNING", "ERROR"]
    parser.add_argument(
        "-l",
        "--log",
        dest="log",
        type=str,
        default=log_levels[1],
        choices=log_levels,
        help=f"The level of log. If it is not specified the default level is {log_levels[1]}"
    )
    parser.add_argument(
        "-g",
        "--gui",
        action='store_true',
        help="Used as one single argument, run the graphical interface application for selecting the arguments and processing the data."
    )
    return parser.parse_args()

@Gooey(program_name='Acceleration profiler application', default_size=(600, 650))
def define_argument_with_gui(default_configuration_file_path: str, signal_profiles: list):
    """Setup the arguments for a graphical interface"""
    parser = GooeyParser(description='This will analyze accelerations')
    required_argument_parser = parser.add_argument_group('Required options')
    required_argument_parser.add_argument(
        "-c",
        "--toml_configuration_file_path",
        required=True,
        dest="toml_configuration_file_path",
        type=str,
        default=default_configuration_file_path
    )
    required_argument_parser.add_argument(
        '--InputMode',
        default='File',
        dest="input_mode",
        choices=['File','Folder'],
        widget='Dropdown'
    )
    required_argument_parser.add_argument(
        '--SignalProfile',
        default="Cadillac_CT6",
        dest="signal_profile",
        choices=signal_profiles,
        widget='Dropdown'
    )
    process_steps = [step.name.lower() for step in list(ProcessStep)]
    required_argument_parser.add_argument(
        '--step',
        default=ProcessStep.PREPROCESSING.name.lower(),
        choices=process_steps,
        dest="process_step",
        help=process_step_help(),
        widget='Dropdown'
    )
    required_argument_parser.add_argument(
        '--MappingFile',
        dest="mapping_file",
        help='File that maps drivers to trips.',
        widget='FileChooser'
    )
    single_group = required_argument_parser.add_argument_group(
        "Single Input Options",
        "Choose if only one file is analyzed"
    )
    single_group.add_argument(
        '--InputFile',
        dest="input_file",
        widget='FileChooser'
    )
    multi_group = required_argument_parser.add_argument_group(
        "Multi Input Options",
        "Choose if an entire folder is analyzed"
    )
    multi_group.add_argument(
        '--InputFolder',
        dest="input_folder_path",
        widget='DirChooser'
    )
    multi_group.add_argument(
        '--UseAll',
        dest="use_all_files",
        help="If True, it applies step on all files in the folder (default: False)",
        action='store_true'
    )
    multi_group.add_argument(
        '--FolderInputStart',
        dest="starting_folder_input",
        default=0,
        widget='IntegerField'
    )
    multi_group.add_argument(
        '--FolderInputEnd',
        default=1,
        dest="ending_folder_input",
        widget='IntegerField'
    )
    parameter_estimation_group = parser.add_argument_group(
        "Parameter Estimation Options",
        "Choose your options for the parameter estimation. These are only relevant if 'modeling' step is chosen."
    )
    estimation_steps = [estimation.value.lower() for estimation in list(EstimationStep)]
    parameter_estimation_group.add_argument(
        '--Estimation_Step',
        dest="estimation_step",
        choices=estimation_steps,
        widget='Dropdown'
    )
    models = [model.value for model in list(Model)]
    parameter_estimation_group.add_argument(
        '--model',
        default=Model.CONSTANT_ACCELERATION_AND_SECOND_ORDER_STEP_RESPONSE.value,
        dest="model",
        choices=models
    )
    matlab_options = parameter_estimation_group.add_argument_group(
        "MATLAB Estimation options",
        "For 'Print_MATLAB_commands'. Requires 'Print_MATLAB_commands'. MATLAB commands automatically saved in clipboard."
    )
    matlab_options.add_argument(
        '--MATLAB_command',
        help="The MATLAB step needs to be specified",
        choices=["parameter_estimation", "save_parameters"],
        dest="matlab_command",
        widget='Dropdown'
    )
    matlab_options.add_argument(
        '--Starting_Speedup_Index',
        dest="starting_speedup_index",
        help='Based on "trip_infos.csv" from "Prepare_MATLAB_Estimation"',
        widget='IntegerField'
    )
    matlab_options.add_argument(
        '--Ending_Speedup_Index',
        dest="ending_speedup_index",
        help='Based on "trip_infos.csv" from "Prepare_MATLAB_Estimation"',
        widget='IntegerField'
    )
    matlab_options.add_argument(
        '--Save_postfix',
        dest="save_postfix",
        help='Postfix with which the file gets saved in the save_parameters command.',
        action='store'
    )
    parameter_optional_group = parser.add_argument_group(
        "Optional parameters",
        "Choose the level of logging"
    )
    log_levels = ["DEBUG", "INFO", "WARNING", "ERROR"]
    parameter_optional_group.add_argument(
        "-l",
        "--log",
        dest="log",
        type=str,
        default=log_levels[1],
        choices=log_levels,
        help=f"The level of log. If it is not specified the default level is {log_levels[1]}"
    )
    return parser.parse_args()

def setup_logging(logging_level: str, trip_name: str, log_directory_path: str):
    """Configure the logging parameters (log output location, level of logging)"""
    current_timestamp = datetime.now().strftime("%Y-%m-%dT%H-%M-%S.%fZ")
    log_file_name = f"{trip_name}-{current_timestamp}.log"
    build_directory_tree([log_directory_path])
    log_file_path = os.path.join(log_directory_path, log_file_name)
    numeric_log_level = getattr(logging, logging_level.upper(), None)
    if not isinstance(numeric_log_level, int):
        raise ValueError('Invalid log level: %s' % logging_level)
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(numeric_log_level)
    stream_handler.setFormatter(logging.Formatter(fmt='%(message)s'))
    logging.basicConfig(
        format='%(asctime)s;%(levelname)s;%(funcName)s;%(message)s',
        level=logging.DEBUG,
        handlers=[
            logging.FileHandler(log_file_path, mode='w'),
            stream_handler
        ]
    )

def main():
    """Main function setting up the logging, the parsing and running the data treatment."""
    default_configuration_file_path = PathConfiguration.default_paths_configuration_file_path()
    with open(signal_profiles_path(), 'r', encoding='utf-8') as signals_file:
        all_signal_profiles = json.load(signals_file)
    profile_name_list = list(all_signal_profiles.keys())[:-1]
    # Use CLI or HMI
    arguments = None
    if '--gui' in sys.argv or '-g' in sys.argv:
        arguments = define_argument_with_gui(
            default_configuration_file_path=default_configuration_file_path,
            signal_profiles=profile_name_list
        )
    else:
        arguments = define_argument_with_cli(
            default_configuration_file_path=default_configuration_file_path,
            signal_profiles=profile_name_list
        )
    arguments.process_step = ProcessStep[arguments.process_step.upper()]
    # Paths
    path_configuration: PathConfiguration = None
    with open(arguments.toml_configuration_file_path, 'rt', encoding='utf-8') as configuration_file:
        path_configuration = PathConfiguration(toml.load(configuration_file))
    # Initiate the configuration
    estimation_settings = {}
    if arguments.process_step.value >= ProcessStep.MODELING.value:
        estimation_settings = {
            "Estimation_Step": arguments.estimation_step,
            "MATLAB_command": arguments.matlab_command,
            "model": arguments.model,
            "start": arguments.starting_speedup_index,
            "end": arguments.ending_speedup_index,
            "postfix": arguments.save_postfix
        }
    configuration = Configuration(
        estimation_settings=estimation_settings,
        subject_trip_mapping_file_path=arguments.mapping_file,
        signal_profile=all_signal_profiles[arguments.signal_profile],
        path_configuration=path_configuration
    )
    trip_files_to_process = []
    print(arguments)
    if arguments.input_mode == "File":
        trip_file = TripFile(file_path=arguments.input_file)
        trip_files_to_process.append(trip_file)
    elif arguments.input_mode == 'Folder':
        folder_path = arguments.input_folder_path
        if arguments.use_all_files:
            for input_file in os.listdir(folder_path):
                trip_file = TripFile(file_path=os.path.join(folder_path, input_file))
                trip_files_to_process.append(trip_file)
        else:
            for i, input_file in enumerate(os.listdir(folder_path)):
                if arguments.starting_folder_input <= i <= arguments.ending_folder_input:
                    trip_file = TripFile(file_path=os.path.join(folder_path, input_file))
                    trip_files_to_process.append(trip_file)
    for trip_file in trip_files_to_process:
        setup_logging(
            logging_level=arguments.log,
            trip_name=trip_file.file_name_without_extension,
            log_directory_path=path_configuration.log_directory_path()
        )
        process_trip_file(
            input_file_path=trip_file.file_path,
            step=arguments.process_step,
            configuration=configuration
        )

if __name__ == "__main__":
    main()
