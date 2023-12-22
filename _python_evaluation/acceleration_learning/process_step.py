"""Python module representing the processing steps"""
from enum import Enum
class ProcessStep(Enum):
    """Represent the acceleration profiler data processing steps"""
    FILE_STRUCTURE_INITIALIZATION = 1 # temp_folder_preparation
    TIME_WINDOW_EXTRACTION = 3  # extraction
    PREPROCESSING = 4
    TITLE_ID_EXTRACTION = 5
    MAP_EXTRACTION = 6
    MAP_ANALYSIS = 7
    BASE_ANALYSIS = 8
    MODELING = 9
