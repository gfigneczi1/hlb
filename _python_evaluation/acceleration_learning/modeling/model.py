"""Python enum to represent the acceleration model"""
from enum import Enum

class Model(Enum):
    """Represent the acceleration model"""
    CONSTANT_ACCELERATION_AND_SECOND_ORDER_STEP_RESPONSE = "CA_PT2"
    SECOND_ORDER_STEP_RESPONSE = "PT2"
    FIRST_ORDER_STEP_RESPONSE = "PT1"
