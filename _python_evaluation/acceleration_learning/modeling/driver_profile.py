from copy import deepcopy
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum

DEFAULT_DICT = {
    DrivingEventEnum.SPEED_UP: 0.5    
}

class DriverProfile:
    def __init__(self, name):
        self.name = name
        self.levels_of_dynamic = deepcopy(DEFAULT_DICT)

    def __str__(self):
        return f"DriverProfile ({self.name}): {self.get_average_dynamic()}"

    def get_average_dynamic(self):
        sum_of_levels = 0
        for level in self.levels_of_dynamic.values():
            sum_of_levels = level
        return sum_of_levels / len(self.levels_of_dynamic.keys())