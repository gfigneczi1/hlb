"""Python module to structure the configuration's variables in one class"""
from dataclasses import dataclass
import toml
import json
from acceleration_learning.path_configuration import PathConfiguration

@dataclass
class Configuration:
    """Configuration parameters container"""
    estimation_settings: dict
    subject_trip_mapping_file_path: str
    signal_profile: dict
    path_configuration: PathConfiguration

def default_configuration() -> Configuration:
    """Set a Configuration with the default paths configuration, the Cadillac signals profile and
    none estimation settings"""
    path_configuration: PathConfiguration = None
    with open(PathConfiguration.default_paths_configuration_file_path(), 'rt', encoding='utf-8') as configuration_file:
        path_configuration = PathConfiguration(toml.load(configuration_file))
    cadillac_signals_profile = {}
    with open(path_configuration.data_signal_profile_file_path(), 'r', encoding='utf-8') as signal_profile_file:
        cadillac_signals_profile = json.load(signal_profile_file)["Cadillac_CT6"]
    return Configuration(
        estimation_settings=None,
        subject_trip_mapping_file_path=path_configuration.driver_trip_mapping_file_path,
        signal_profile=cadillac_signals_profile,
        path_configuration=path_configuration,
    )
