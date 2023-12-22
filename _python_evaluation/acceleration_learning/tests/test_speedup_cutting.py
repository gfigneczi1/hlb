
import unittest
from acceleration_learning.tests.situation_extracted_period_tester import SituationExtractedPeriodTester
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.preprocessing.situation_extraction.speed_up_cutter import extract_speed_ups

class TestSpeedupCutter(SituationExtractedPeriodTester):

    def test_speed_up_detection(self):
        """Test: verify that the expected time interval of speed ups from reference test file are included in the the extracted speed ups"""
        self.check_detection(
            tested_driving_situation=DrivingEventEnum.SPEED_UP,
            situation_extraction_function=extract_speed_ups,
            ground_truth_file_configuration_name='speedup_ground_truth_file_path'
        )

if __name__ == '__main__':
    unittest.main()
