import pandas as pd
import numpy as np
import os
import sys
import json
import pathlib
import logging
import pyperclip
from scipy.io import loadmat
from acceleration_learning.subject_identifier import SubjectIdentifier

# Define CA-PT2 parameter (and PT2) filter thresholds
D_LOWER_THRESHOLD = 0.25
D_HIGHER_THRESHOLD = 2.0
W0_LOWER_THRESHOLD = 0.0
W0_HIGHER_THRESHOLD = 1.0
TC_LOWER_THRESHOLD = -2.5
TC_HIGHER_THRESHOLD = 10.0
ALFA_LOWER_THRESHOLD = 0.0
ALFA_HIGHER_THRESHOLD = 0.4
# Define PT1 parameter
TP1_LOWER_THRESHOLD = 0.0
TP1_HIGHER_THRESHOLD = 20.0

class matlab_estimation_helper():

    def prepare_matlab_estimation(
        self,
        mapping_file_path: str,
        input_multiple_trip_windows_directory_path: str,
        input_output_postfixed_trip_information_file_path: str
    ):
        '''creates a csv file called "trip_infos_postfix.csv" 
        that includes meta information mapping the trip to its driver_id and the number of speedups in it.
        
        Arguments:
        input_preprocessed_windows_directory_path -- the path to the directory containing all the trip directories
        '''
        data_folder_path = input_multiple_trip_windows_directory_path
        if os.path.exists(input_output_postfixed_trip_information_file_path):
            df = pd.read_csv(input_output_postfixed_trip_information_file_path)
            names = list(df["trip_name"])
            number_of_speed_ups_in_folders = list(df["number_of_speedups"])
            subjects = list(df["subject_id"])
        else:
            names = []
            number_of_speed_ups_in_folders = []
            subjects = []
        for folder in os.listdir(data_folder_path):
            number_of_speed_ups = 0
            if not "." in folder and folder not in names:
                for file in os.listdir(os.path.join(data_folder_path, folder)):
                    if "speed" in file:
                        number_of_speed_ups += 1
                if number_of_speed_ups > 0:
                    logging.info("added: %s", folder)
                    names.append(folder)
                    number_of_speed_ups_in_folders.append(number_of_speed_ups)
                    subjects.append(str(SubjectIdentifier.get_subject_id(
                        trip_name=folder,
                        mapping_file_path=mapping_file_path
                    )))
        all_infos = np.array([names, number_of_speed_ups_in_folders, subjects]).T
        df = pd.DataFrame(all_infos, columns=["trip_name", "number_of_speedups", "subject_id"])
        df.to_csv(input_output_postfixed_trip_information_file_path)
        logging.info("The MATLAB estimation can be done now!")

    def print_matlab_commands(
        self,
        estimation_settings,
        signal_profile,
        postfixed_matlab_model_parameters_file_path: str,
        input_postfixed_trip_information_file_path: str,
        extraction_output_directory_path: str
        ):
        '''prints out the desired MATLAB commands and copies them into the clipboard.
        This happens if "modeling" is chosen in "step" and "Print_MATLAB_commands" is chosen in "Estimation step"

        keyword arguments:
        estimation_settings -- settings set in the GUI in a json files with keys: 
                                    - "Estimation_Step" (not relevant in this method)
                                    - "Python_Optimizer" (not relevant in this method)
                                    - "MATLAB_command" ("parameter_estimation" or "save_parameters")
                                    - "model" (CA_PT2, PT2, PT1)
                                    - "start" (start in trip_infos_postfix.csv)
                                    - "end" (end in trip_infos_postfix.csv)
                                    - "postfix" (postfix attached to file that saves the estimated parameters)    
        '''
        trip_info_path = input_postfixed_trip_information_file_path
        copy_paste_text = ""
        if not estimation_settings["model"] in ["PT2", "PT1","CA_PT2"]:
            sys.exit()
        if os.path.exists(trip_info_path):
            df = pd.read_csv(trip_info_path)
            if estimation_settings["MATLAB_command"] == "parameter_estimation":
                logging.info("Matlab instructions to paste into the MATLAB console for estimating parameters:")
                end = min(df.shape[0], estimation_settings["end"])
                for i in range(estimation_settings["start"], end):
                    file_name = df["trip_name"][i]
                    file_number = df["number_of_speedups"][i]
                    for j in range(file_number):
                        model = estimation_settings["model"]
                        time = signal_profile["Time"]
                        longitudinal_velocity = signal_profile["Longitudinal_Velocity"]
                        trip_extraction_output_directory_path = os.path.join(extraction_output_directory_path, file_name)
                        matlab_command = f"trip{i}_speedup_{j}_parameter = individual_speedup_parameter_estimation(\"{trip_extraction_output_directory_path}\", {j}, \"{model}\", \"{time}\", \"{longitudinal_velocity}\");"
                        logging.info(matlab_command)
                        copy_paste_text += matlab_command
            elif estimation_settings["MATLAB_command"] == "save_parameters":
                copy_paste_text = f"save(\"{postfixed_matlab_model_parameters_file_path}\")"
                logging.info("Matlab instruction to paste into the MATLAB console for saving the MATLAB environment:")
                logging.info(copy_paste_text)
            else:
                logging.info("what do you want to do with MATLAB? Specify the mode from \n 'parameter_estimation'  \n 'save_parameters'")
            pyperclip.copy(copy_paste_text)
    
    def get_trip_info(self, window_information_file_path: str) -> list:
        """Return the list of time window details for the given trip file """
        with open(window_information_file_path, 'r') as open_file:
            read_list = json.load(open_file)
        return read_list
    
    def write_trip_info(self, write_list: list, window_information_file_path):
        """Save the window trip information list to the given path
        (overwrite existing file)"""
        with open(window_information_file_path, 'w', encoding="utf-8") as new_file:
            json.dump(write_list, new_file, indent=4)

    def filter_out_parameters(self, parameters, model):
        '''this method filters out parameters that are statistically outliers'''
        if model == "PT1":
            if parameters[1] < TP1_LOWER_THRESHOLD or parameters[1] > TP1_HIGHER_THRESHOLD:
                return True
            else:
                return False
        if model == "PT2":
            if parameters[1] < W0_LOWER_THRESHOLD or parameters[1] > W0_HIGHER_THRESHOLD:
                return True
            elif parameters[2] < D_LOWER_THRESHOLD or parameters[2] > D_HIGHER_THRESHOLD:
                return True
            else:
                return False
        else:
            if parameters[0] < D_LOWER_THRESHOLD or parameters[0] > D_HIGHER_THRESHOLD:
                return True
            elif parameters[1] < W0_LOWER_THRESHOLD or parameters[1] > W0_HIGHER_THRESHOLD:
                return True
            elif parameters[2] < TC_LOWER_THRESHOLD or parameters[2] > TC_HIGHER_THRESHOLD:
                return True
            elif parameters[3] < ALFA_LOWER_THRESHOLD or parameters[3] > ALFA_HIGHER_THRESHOLD:
                return True
            else:
                return False

    def save_trip_parameter_output(self, trip_parameters, model, trip_parameter_output_file_path: str):
        """Save the given model parameters to the given file as a JSON file
        (overwrite existing file)"""
        if model == "PT1":
            parameter_dictionary = {"K": trip_parameters[0], "Tp1": trip_parameters[1]}
        elif model == "PT2":
            parameter_dictionary = {"K": trip_parameters[0], "D": trip_parameters[1], "w0": trip_parameters[2]}
        else:
            parameter_dictionary = {
                "K": 1.0,
                "D": trip_parameters[0],
                "w0": trip_parameters[1],
                "tc": trip_parameters[2],
                "alfa": trip_parameters[3],
                "t_delay": trip_parameters[4]
            }
        with open(trip_parameter_output_file_path, 'w', encoding="utf-8") as new_file:
            json.dump(parameter_dictionary, new_file, indent=4)

    def save_matlab_estimates_to_info_files(
        self,
        estimation_settings,
        postfixed_matlab_model_parameters_file_path: str,
        postfixed_trip_information_file_path: str,
        window_information_file_path: str,
        trip_parameter_output_file_path: str,
        ):
        '''This saves the matlab estimates to the metadata info json files'''        
        mat_dict = loadmat(postfixed_matlab_model_parameters_file_path)
        df = pd.read_csv(postfixed_trip_information_file_path)
        
        end = min(df.shape[0], estimation_settings["end"])
        all_trip_parameters = []
        logging.info("Saving all parameters predicted in MATLAB.")
        for i in range(estimation_settings["start"], end):
            file_name = df["trip_name"][i]
            file_number = df["number_of_speedups"][i] - 1
            read_list = self.get_trip_info(
                window_information_file_path=window_information_file_path
            )
            write_list = []
            j = 0
            for read_window_dict in read_list:
                if read_window_dict["classification"] == "speed_up":
                    parameter_string = "trip{index}_speedup_{speedup_index}_parameter".format(index = i, speedup_index = j)
                    parameters = mat_dict[parameter_string][0]
                    write_window_dict = read_window_dict.copy()
                    write_window_dict["modeling_parameters"] = list(parameters)
                    if not self.filter_out_parameters(parameters, estimation_settings["model"]):
                        all_trip_parameters.append(list(parameters))
                    j += 1
                else:
                    write_window_dict = read_window_dict.copy()
                write_list.append(write_window_dict)
            self.write_trip_info(
                write_list=write_list,
                window_information_file_path=window_information_file_path
            )
            self.save_trip_parameter_output(
                trip_parameters=np.array(all_trip_parameters).mean(axis=0),
                model=estimation_settings["model"],
                trip_parameter_output_file_path=trip_parameter_output_file_path
            )
        logging.info("Saving ended.")
