import unittest
from acceleration_learning.tests.situation_extracted_period_tester import SituationExtractedPeriodTester
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.preprocessing.situation_extraction.slow_down_cutter import extract_slow_downs

class TestSlowdownCutter(SituationExtractedPeriodTester):

    def test_slow_down_detection(self):
        """Test: verify that the expected time interval of slow downs from reference test file are included in the the extracted slow downs"""
        self.check_detection(
            tested_driving_situation=DrivingEventEnum.SLOW_DOWN,
            situation_extraction_function=extract_slow_downs,
            ground_truth_file_configuration_name='slowdown_ground_truth_file_path'
        )

if __name__ == '__main__':
    unittest.main()
