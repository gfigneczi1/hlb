"""Python module to help testing driving situation detection"""
import os
import pathlib
from typing import List, Callable, Tuple
import toml
import json
import unittest
import pandas as pd
from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.path_configuration import PathConfiguration

class SituationExtractedPeriodTester(unittest.TestCase):
    
    def setUp(self):
        with open(PathConfiguration.default_paths_configuration_file_path(), 'rt', encoding='utf-8') as config_file:
            self.path_configuration = PathConfiguration(toml.load(config_file))
        self.data_root_directory = self.path_configuration.preprocessing_input_directory()
    
    def check_detection(self, tested_driving_situation: DrivingEventEnum, situation_extraction_function:Callable[[pd.DataFrame, str, dict], Tuple[list, list, list]], ground_truth_file_configuration_name: str):
        """Test: verify that the expected time interval of the given situation from reference test file are included in the the extracted time windows"""
        with open(self.path_configuration.data_signal_profile_file_path(), 'rt', encoding='utf-8') as profile_file:
            signal_profile = json.load(profile_file)['Cadillac_CT6']
        expected_timing_file_path: str = self.path_configuration.path_configuration[ground_truth_file_configuration_name]
        expected_situation: DrivingEventEnum = tested_driving_situation
        with open(expected_timing_file_path, encoding='utf-8') as expected_timing_file_json:
            expected_timings = json.load(expected_timing_file_json)
            total_expected_driving_situation_n = 0
            total_detected_driving_situation_n = 0
            print()
            for trip_properties in expected_timings.values():
                # For one trip
                filename: str = trip_properties["filename"]
                filepath: str = os.path.join(self.data_root_directory, filename)
                data_frame: pd.DataFrame = pd.read_csv(filepath)
                # Cut situation windows
                driving_situations: [driving_situation_window]
                driving_situations, _, _ = situation_extraction_function(data_frame, filename, signal_profile)
                n_expected_driving_situations = len(trip_properties["time"])
                total_expected_driving_situation_n += n_expected_driving_situations
                n_detected_driving_situations = 0
                for expected_start_and_end in trip_properties["time"]:
                    # For one interval
                    with self.subTest(f"In the file {filename}, for {expected_situation} situation, expected time interval: {expected_start_and_end}"):
                        current_situation_detected: bool = False
                        tested_windows_log: list(str) = []
                        for tested_window in driving_situations:
                            tested_windows_log.append(tested_window)
                            tested_time_interval = tested_window.get_time_interval()
                            if expected_start_and_end[0] < tested_time_interval[0]:
                                break
                            if tested_time_interval[0] <= expected_start_and_end[0] and expected_start_and_end[1] <= tested_time_interval[1]:
                                current_situation_detected = True
                                break
                        if current_situation_detected:
                            n_detected_driving_situations += 1
                            self.assertTrue(True)
                        else:
                            error_message = f"Not detected: expected {expected_situation.value} {str(expected_start_and_end)} in the file {filepath}"
                            n_last = 2
                            if len(tested_windows_log) > n_last:
                                error_message += f"\n the {n_last} last closest window(s) was/were:"
                                for i in range(len(tested_windows_log)-1, max(0, len(tested_windows_log)-1-n_last), -1):
                                    if n_last == 0:
                                        break
                                    error_message += f"\n log: {tested_windows_log[i]}"
                                    n_last -= 1
                            self.assertTrue(False, msg=error_message)
                total_detected_driving_situation_n += n_detected_driving_situations
                if n_expected_driving_situations:
                    detection_ratio = n_detected_driving_situations / n_expected_driving_situations
                    print(f"In the file {filename}, {expected_situation.value} detection ratio is {detection_ratio:.2f} ({n_detected_driving_situations}/{n_expected_driving_situations})")
        total_ratio = total_detected_driving_situation_n / total_expected_driving_situation_n
        print(f"Total detection score: {total_detected_driving_situation_n}/{total_expected_driving_situation_n} ({total_ratio:.2f})")