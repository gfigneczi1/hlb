"""Python class to represent a trip file"""
import os

class TripFile:
    """Simple file representation with system path and name"""
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.file_name = file_path.split(str(os.sep))[-1]
        self.file_name_without_extension = self.file_name.split(os.extsep)[0]
