"""Python unittest module to check the TimeWindowExtractor class"""
import os
import pathlib
import unittest
import json
import toml
import pandas as pd

from acceleration_learning.preprocessing.situation_extraction.time_window_extractor import TimeWindowExtractor
from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.path_configuration import PathConfiguration


def tested_bounds_include_expected_bounds(tested_bounds, expected_bounds):
    """ Return True if the expected_bounds values include the tested_bounds's ones with a tolerance
    
    Keyword arguments:
    expected_bounds -- list of the start and the end index of the expected values (2 elements expected)
    tested_bounds -- list of the start and the end index of the tested values (2 elements expected)
    """
    TIME_ACCURACY_TOLERANCE_IN_SEC = 3
    return expected_bounds[0] - TIME_ACCURACY_TOLERANCE_IN_SEC <= tested_bounds[0] and\
        tested_bounds[0] <= expected_bounds[0] and\
        expected_bounds[1] <= tested_bounds[1] and\
        tested_bounds[1] <= expected_bounds[1] + TIME_ACCURACY_TOLERANCE_IN_SEC

class TestWindowMaker(unittest.TestCase):
    def setUp(self):
        with open(PathConfiguration.default_paths_configuration_file_path(), 'rt', encoding='utf-8') as config_file:
            self.path_configuration = PathConfiguration(toml.load(config_file))
        self.data_root_directory = self.path_configuration.preprocessing_input_directory()
    
    def assist_test_driving_situation_detection(self, situation: DrivingEventEnum, ground_truth_file_configuration_name: str):
        """Check the timing of the extracted windows of only the given situation"""
        expected_timing_file_path: str = self.path_configuration.path_configuration[ground_truth_file_configuration_name]
        data_root_file: str = self.data_root_directory
        with open(self.path_configuration.data_signal_profile_file_path(), 'rt', encoding='utf-8') as profile_file:
            cadillac_signal_profile = json.load(profile_file)['Cadillac_CT6']
        with open(expected_timing_file_path, 'rt', encoding='utf-8') as driving_situation_file:
            any_situation_per_trip = json.load(driving_situation_file)
            for trip_value in any_situation_per_trip.values():
                if situation.value not in trip_value:
                    continue
                filepath = os.path.join(data_root_file, trip_value["filename"])
                data_frame = pd.read_csv(filepath)
                situation_windows = TimeWindowExtractor(trip_name=filepath).makeWindows(data_frame, cadillac_signal_profile)
                extracted_situations = situation_windows[situation]
                for expected_start_and_end in trip_value[situation.value]:
                    current_situation_detected = False
                    for tested_window in extracted_situations:
                        tested_time_interval = [tested_window.data[cadillac_signal_profile["Time"]].iloc[0], tested_window.data[cadillac_signal_profile["Time"]].iloc[-1]]
                        if tested_bounds_include_expected_bounds(tested_time_interval, expected_start_and_end):
                            current_situation_detected = True
                            break
                    if current_situation_detected:
                        self.assertTrue(True)
                    else:
                        self.assertTrue(False, msg=f"The expected {situation.value} {str(expected_start_and_end)} in the file {filepath}  is not detected")

    def test_all_situations_detection(self):
        """Test: verify every time interval for every situations in the reference test file are detected by the windows maker"""
        expected_timing_file_path: str = self.path_configuration.path_configuration['any_situation_ground_truth_file_path']
        data_root_file: str = self.data_root_directory
        with open(self.path_configuration.data_signal_profile_file_path(), 'rt', encoding='utf-8') as profile_file:
            cadillac_signal_profile = json.load(profile_file)['Cadillac_CT6']
        with open(expected_timing_file_path, 'rt', encoding='utf-8') as expected_timing_file_json:
            expected_timings = json.load(expected_timing_file_json)
            detection_ratio_per_situation = {}
            for situation in DrivingEventEnum:
                detection_ratio_per_situation[situation] = {
                    "n_detected_situations": 0,
                    "n_expected_situations": 0,
                }
            print()
            for trip_properties in expected_timings.values():
                filename: str = trip_properties["filename"]
                filepath: str = os.path.join(data_root_file, filename)
                data_frame: pd.DataFrame = pd.read_csv(filepath)
                wm = TimeWindowExtractor(trip_name=filename)
                situation_windows: dict = wm.makeWindows(
                    df=data_frame,
                    signal_profile=cadillac_signal_profile
                )
                expected_situations: dict = trip_properties["situations"]
                for expected_situation, expected_time_intervals in expected_situations.items():
                    expected_situation: DrivingEventEnum = DrivingEventEnum(expected_situation)
                    filtered_speed_ups: list = situation_windows[expected_situation]
                    filtered_speed_ups.sort()
                    n_expected_situations = len(expected_time_intervals)
                    n_detected_situations = 0
                    for expected_start_and_end in expected_time_intervals:
                        with self.subTest(f"In file {filename}, for {expected_situation} situation, expected time interval: {expected_start_and_end}"):
                            current_situation_detected: bool = False
                            tested_windows_log = []
                            for tested_window in filtered_speed_ups:
                                tested_windows_log.append(tested_window)
                                tested_time_interval = tested_window.get_time_interval()
                                if expected_start_and_end[0] < tested_time_interval[0]:
                                    break
                                if tested_bounds_include_expected_bounds(tested_time_interval, expected_start_and_end):
                                    current_situation_detected = True
                                    break
                            if current_situation_detected:
                                self.assertTrue(True)
                                n_detected_situations += 1
                            else:
                                error_message = f"Not detected: expected {expected_situation.value} {str(expected_start_and_end)} in the file {filepath}"
                                n_last = 1
                                if len(tested_windows_log) > n_last:
                                    error_message += f"\n the {n_last} last closest windows were:"
                                    for i in range(len(tested_windows_log)-1, 0, -1):
                                        if n_last == 0:
                                            break
                                        error_message += f"\n log: {tested_windows_log[i]}"
                                        n_last -= 1
                                self.assertTrue(False, msg=error_message)
                    if n_expected_situations:
                        detection_ratio = n_detected_situations / n_expected_situations
                        print(f"In the file {filename}, {expected_situation} detection ratio is {detection_ratio:.2f} ({n_detected_situations}/{n_expected_situations})")
                    detection_ratio_per_situation[expected_situation]['n_detected_situations'] += n_detected_situations
                    detection_ratio_per_situation[expected_situation]['n_expected_situations'] += n_expected_situations
            print()
            for situation, detected_expected in detection_ratio_per_situation.items():
                if detected_expected['n_expected_situations'] == 0:
                    ratio = 0
                else:
                    ratio = detected_expected['n_detected_situations'] / detected_expected['n_expected_situations']
                print(f"Total detection score of {situation.value}: {detected_expected['n_detected_situations']}/{detected_expected['n_expected_situations']} ({ratio:.2f})")
            

    @unittest.skip("Skipped because this test is included in the all_situations_detection test")
    def test_situations(self):
        """Test: check if the identified change lanes in the reference JSON file are detected by the windows maker"""
        self.assist_test_driving_situation_detection(DrivingEventEnum.LANE_CHANGE)

    @unittest.skip("Skipped because this test is included in the all_situations_detection test")
    def test_speed_up_detection(self):
        """Test: check if the identified speed up in the JSON file are detected by the windows maker"""
        self.assist_test_driving_situation_detection(DrivingEventEnum.SPEED_UP)

    @unittest.skip("Skipped because this test is included in the all_situations_detection test")
    def test_slow_down_detection(self):
        """Test: check if the identified slow down in the JSON file are detected by the windows maker"""
        self.assist_test_driving_situation_detection(DrivingEventEnum.SLOW_DOWN)

if __name__ == '__main__':
    unittest.main()
