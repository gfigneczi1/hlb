"""Python module to represent the project file structure configuration"""
import os
from acceleration_learning.trip_file import TripFile
from utils.path import configuration_directory_path, root_path

class PathConfiguration():
    """Class in charge of storing the directory locations for input, temporary and output files"""
    WORKING_DIRECTORY_KEYWORD = 'working_directory'
    WORKING_DIRECTORY_RELATIVE_PATH_KEYWORD = 'path_is_relative_to_the_working_directory'
    PROJECT_ROOT_RELATIVE_PATH_KEYWORD = 'path_is_relative_to_the_project_root_directory'
    KEYWORDS = [
        WORKING_DIRECTORY_RELATIVE_PATH_KEYWORD,
        PROJECT_ROOT_RELATIVE_PATH_KEYWORD
    ]

    @classmethod
    def default_paths_configuration_file_path(cls) -> str:
        """Return the default path to the configuration file"""
        return os.path.join(configuration_directory_path(), "default_paths_configuration.toml")

    def __init__(self, toml_configuration_dict: dict):
        configuration = dict()
        working_directory_path = toml_configuration_dict[self.WORKING_DIRECTORY_KEYWORD]
        for toml_table, toml_pairs in toml_configuration_dict.items():
            path_prefix = ""
            if toml_table == self.WORKING_DIRECTORY_KEYWORD:
                continue
            if toml_pairs.get(self.WORKING_DIRECTORY_RELATIVE_PATH_KEYWORD):
                path_prefix = working_directory_path
            elif toml_pairs.get(self.PROJECT_ROOT_RELATIVE_PATH_KEYWORD):
                path_prefix = root_path()
            for path_keyword, path_value in toml_pairs.items():
                if path_keyword in self.KEYWORDS:
                    continue
                configuration[path_keyword] = os.path.join(path_prefix, path_value)
        self.path_configuration = configuration

    def directory_hierarchy_string(self) -> str:
        """Return the project file hierarchy """
        paths = list(self.path_configuration.values())
        paths.sort()
        return "\n".join(paths)

    def log_directory_path(self) -> str:
        """Return the path to the log directory"""
        return self.path_configuration["log_directory"]

    def data_signal_profile_file_path(self) -> str:
        """Return the path to the signal profile file"""
        return self.path_configuration["signal_profile_file"]
    
    def driver_trip_mapping_file_path(self) -> str:
        """Return the path to the driver to trip mapping file"""
        return self.path_configuration['subject_trip_mapping_file']
    
    def preprocessing_input_directory(self) -> str:
        """Return the path to the trip data directory"""
        return self.path_configuration["preprocessing_input_directory"]

    def preprocessing_extraction_output_directory_path(self) -> str:
        """Return the output directory path of the window time extraction
        step"""
        return self.path_configuration["preprocessing_extraction_output_directory"]

    def trip_file_extraction_output_directory_path(self, trip_file: TripFile) -> str:
        """Return the path to the directory storing the extraction's output"""
        return os.path.join(
            self.preprocessing_extraction_output_directory_path(),
            trip_file.file_name_without_extension
        )

    def preprocessing_extraction_information_output_directory_path(self) -> str:
        """Return the path to the directory storing meta information of the 
        extracted time windows"""
        return self.path_configuration["preprocessing_extraction_information_output_directory"]

    def preprocessing_extraction_information_output_file_path(self, trip_file: TripFile) -> str:
        """Return the path to the directory storing meta information of the
        extracted time windows"""
        return os.path.join(
            self.preprocessing_extraction_information_output_directory_path(),
            f"{trip_file.file_name_without_extension}_info.json"
        )

    def preprocessing_map_output_directory_path(self) -> str:
        """Return the path to the directory to store generated map files"""
        return self.path_configuration["preprocessing_map_output_directory"]

    def geographic_data_output_directory_path(self) -> str:
        """Return the path to the output directory of the GeoJson files"""
        return self.path_configuration["geographic_data_output_directory"]

    def geographic_trip_directory_output_directory_path(self, trip_file: TripFile) -> str:
        """Return the path to the directory storing the windows including map information"""
        return os.path.join(
            self.geographic_data_output_directory_path(),
            trip_file.file_name_without_extension
        )

    def geographic_data_output_file_path(self, trip_file: TripFile) -> str:
        """Return the path to the output file of the GeoJson files"""
        return os.path.join(
            self.geographic_data_output_directory_path,
            f"{trip_file.file_name_without_extension}.geojson"
        )

    def map_info_path(self) -> str:
        """Return the path to the map info file"""
        return os.path.join(
            self.preprocessing_map_output_directory_path(),
            "map_info.csv"
        )

    def modeling_output_directory_path(self) -> str:
        """Return the directory storing the model parameters"""
        return self.path_configuration["modeling_output_directory"]

    def modeling_output_trip_parameter_path(self, trip_file: TripFile) -> str:
        """Return the path to the file to store the trip parameters"""
        return os.path.join(
            self.modeling_output_directory_path(),
            f"{trip_file.file_name_without_extension}_parameters.json"
        )

    def modeling_matlab_directory_path(self) -> str:
        """Return the path to the directory to store matlab necessary files"""
        return self.path_configuration["modeling_matlab_directory"]

    def modeling_matlab_model_parameters_input_file_path(self, postfix_name: str) -> str:
        """Return the path to the MATLAB file containing models parameters"""
        return os.path.join(
            self.modeling_matlab_directory_path(),
            f"_all_parameters_{postfix_name}.mat"
        )

    def modeling_trip_information_input_file_path(self, postfix_name: str) -> str:
        """Return the path to a CSV file mapping trip name, subject and number of speed ups"""
        return os.path.join(
            self.modeling_matlab_directory_path(),
            f"trip_infos_{postfix_name}.csv"
        )

    def analysis_weighted_formula_output_directory_path(self) -> str:
        """Return the path to the directory to store weighted formula intermediate and results
        outputs"""
        return os.path.join(
            self.path_configuration["analysis_weighted_formula_output_directory"]
        )
        
    def analysis_curvature_output_directory_path(self) -> str:
        """Return the path to the directory storing the results of the curvature analysis"""
        return os.path.join(
            self.path_configuration["analysis_curvature_output_directory"]
        )
