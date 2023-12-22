import os
import json
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from utils.path import root_path
from acceleration_learning.preprocessing.filtering.filter import filter_acc, filter_non_free_driving
from acceleration_learning.preprocessing.situation_extraction.speed_change_tools import merge_close_windows

def setup_argument_parser() -> argparse.ArgumentParser:
    """Return the configured arguments parser"""
    DEFAULT_LOG_FOLDER: str = os.path.join(root_path(), "_temp", "log")
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="Parse the application log to make statistics of the speed ups extraction"
    )
    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument(
        "-n",
        "--name",
        required=True,
        dest="filename",
        type=str,
        help="The name of the log file."
    )
    parser.add_argument(
        "-d",
        "--directory",
        dest="log_directory_path",
        type=str,
        default=DEFAULT_LOG_FOLDER,
        help="The path to the log directory."
    )
    return parser

def parse_log_file(file_path: str):
    line_separator_character = ";"
    merged_windows_rows = []
    acc_filtered_windows_rows = []
    non_free_driving_windows_rows = []
    with open(file_path, "rt", encoding="utf-8") as log_file:
        for line in log_file:
            line_without_carrier_return = line.rstrip("\n")
            log_line_fields = line_without_carrier_return.split(line_separator_character)
            log_level = log_line_fields[1]
            if log_level != "DEBUG":
                continue
            _, log_level, function_name, log_message = log_line_fields
            log_dict = json.loads(log_message)
            if function_name == merge_close_windows.__name__:
                row = log_dict["merged_windows"]
                interval = {
                    "trip_name": log_dict["name"],
                    "gap_duration_in_sec": row[1]["time_interval"][0] - row[0]["time_interval"][1]
                }
                merged_windows_rows.append(interval)
            elif function_name == filter_acc.__name__:
                interval = {
                    "trip_name": log_dict["name"],
                    "window_time_duration_in_sec": log_dict["time_interval"][1] - log_dict["time_interval"][0]
                }
                acc_filtered_windows_rows.append(interval)
            elif function_name == filter_non_free_driving.__name__:
                interval = {
                    "trip_name": log_dict["name"],
                    "window_time_duration_in_sec": log_dict["time_interval"][1] - log_dict["time_interval"][0],
                    "non_free_driving_duration_in_sec": log_dict["non_free_driving_duration_in_s"]
                }
                non_free_driving_windows_rows.append(interval)
    merged_windows_df = pd.DataFrame(merged_windows_rows)
    acc_filtered_windows_df = pd.DataFrame(acc_filtered_windows_rows)
    non_free_driving_windows_df = pd.DataFrame(non_free_driving_windows_rows)
    return {
        "merged_windows_df": merged_windows_df,
        "acc_filtered_windows_df": acc_filtered_windows_df,
        "non_free_driving_windows_df": non_free_driving_windows_df
    }

def main(args) -> None:
    file_path = os.path.join(args.log_directory_path, args.filename)
    dfs = parse_log_file(file_path=file_path)

    for parsed_message_type, rows in dfs.items():
        database_for_current_parsed_message_file_name = f"{parsed_message_type}.csv"
        database_for_current_parsed_message_file_path = os.path.join(
            root_path(),
            database_for_current_parsed_message_file_name
        )
        file = open(database_for_current_parsed_message_file_path, "w")
        file.close()
        database_df = pd.read_csv(database_for_current_parsed_message_file_path)
        database = database_df.to_dict("index")
        columns: list(str) = []
        first_value = database.values[0]
        if first_value:
            columns = list(first_value.keys())
        for i_row, new_row in enumerate(rows):
            new_dict_row = {}
            for column in columns:
                new_dict_row[column] = new_row[column]
            database[i_row] = new_dict_row
        updated_database_df = pd.DataFrame(database)
        updated_database_df.to_csv(database_for_current_parsed_message_file_path)

if __name__ == "__main__":
    parser: argparse.ArgumentParser = setup_argument_parser()
    main(parser.parse_args())
