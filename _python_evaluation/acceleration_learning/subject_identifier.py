import requests
import os
import pandas as pd
from utils.path import root_path

DEFAULT_SUBJECT_FILE_PATH = f"{root_path()}\\_temp\\subjects.csv"
class SubjectIdentifier:
    """Associate a trip name to its driver"""
    def __init__(self, subject_file_path=DEFAULT_SUBJECT_FILE_PATH):
        self.subject_file_path = subject_file_path
        self.subject_df = self.load_subject_file(subject_file_path)

    def get_driver_id(self, trip_label: str) -> str:
        """Return the subject number corresponding to the given trip"""
        return self.subject_df[self.subject_df["trip_label"] == trip_label]["subject_id"].iloc[0]

    def load_subject_file(self, download_destination_file_path: str) -> pd.DataFrame:
        """Return the local CSV file as a data frame"""
        return pd.read_csv(download_destination_file_path, sep=";")

    @staticmethod
    def get_subject_id(trip_name: str, mapping_file_path: str) -> str:
        '''returns the id of a subject out of the Mapping file chosen in the GUI
        Make sure to download an up-to-date version of this mapping file from the FDAL repository

        keyword arguments:
        trip_name -- name of the trip that is analyzed
        '''
        if os.path.exists(mapping_file_path):
            data = pd.read_csv(mapping_file_path, sep=";")
            for trip_label, subject_id in zip(data["trip_label"], data["subject_id"]):
                if trip_label in trip_name:
                    return subject_id
        return "unknown"
