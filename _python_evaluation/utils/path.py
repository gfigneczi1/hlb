"""Help module to generate path to project directories"""
import subprocess
import os
from typing import List

def root_path() -> str:
    """Return the path to the git project root"""
    result = subprocess.run(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE)
    path = result.stdout.decode("utf-8")
    path = path.replace("\n", "")
    path = path.replace('/', '\\')
    return path

def temp_path() -> str:
    """Return the path to the temp folder of the git project"""
    return os.path.join(root_path(), "_temp")

def configuration_directory_path() -> str:
    """Return the path to the configuration folder of the PSN acceleration project"""
    return os.path.join(root_path(), "_python_evaluation", "configuration")

def signal_profiles_path() -> str:
    """Return the path to the JSON file containing the signal profiles of different data format"""
    return os.path.join(
        configuration_directory_path(),
        "acceleration_learning_signal_profiles.json"
    )

def get_csv_file_path(data_folder_path: str, file_index: str):
    ''' returns the path of the trip file that is analyzed (filePath) and just
        the name of the file that is analyzed without the extension (filename)    
    '''
    file_name = os.listdir(data_folder_path)[file_index]
    file_path = os.path.join(data_folder_path, file_name)
    return file_path, file_name

def build_directory_tree(paths: List[str]):
    """Create the paths in the project directory"""
    for path in paths:
        if not os.path.exists(path) and not os.path.isfile(path):
            os.makedirs(path)
