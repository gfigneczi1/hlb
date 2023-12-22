import os
import pathlib
from pathlib import Path
from _python_evaluation.utils.utilities import utilities

class file_getter(): 
    '''
    All of the files in the _temp folder that have the .mat extension are automatically added to the file_list and moved to 0_data_mat
    '''
    def get_automatic_list(self):
        file_list = []
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        top = str(pathlib.Path(__file__).parent.parent.parent)
        paths = [os.path.join(py_eval_path, '0_data_mat')]
        path = os.path.join(top, '_temp')
        util = utilities()
        util.create_path(paths)
        for file in os.listdir(path):
            if '.mat' in file:
                file_list.append(file[:-4])
                os.replace(os.path.join(path, file), os.path.join(py_eval_path, '0_data_mat', file))
        return file_list