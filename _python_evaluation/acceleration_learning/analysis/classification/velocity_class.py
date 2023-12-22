"""Module to represent the speed up velocity classes"""
from enum import Enum

class VelocityClass(Enum):
    """The categories to represent speed change depending on the velocity at
    the start and the end of speed up"""
    ZERO_TO_LOW = "zero_to_low_speed_up"
    ZERO_TO_MID = "low_to_mid_speed_up"
    ZERO_TO_HIGH = "zero_to_high_speed_up"
    LOW_TO_LOW = "low_to_low_speed_up"
    LOW_TO_MID = "low_to_mid_speed_up"
    LOW_TO_HIGH = "low_to_high_speed_up"
    MID_TO_MID = "mid_to_mid_speed_up"
    MID_TO_HIGH = "mid_to_high_speed_up"
    HIGH_TO_HIGH = "high_to_high_speed_up"
    UNKNOWN = "unknown"
