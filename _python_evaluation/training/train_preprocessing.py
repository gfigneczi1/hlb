import scipy.io
import pandas as pd
import os
import glob
import json
from matplotlib import pyplot as plt
import numpy as np
import pathlib
from _python_evaluation.utils.Convert_LaneChange import conversion
from _python_evaluation.utils.Extract_LaneChange import cutter
from _python_evaluation.utils.Relabel_LaneChange import relabeler
from _python_evaluation.utils.file_getter import file_getter

class Train_Preprocessing:
    def execute():
        print("start preprocessing")
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        remove_list = ['__header__', '__version__', '__globals__']
        fg = file_getter()
        file_list = fg.get_automatic_list()

        # not_in_the_file: SteeringTorque, KF_status
        fObj = open(py_eval_path+'/signal_profiles.json')
        signals = json.load(fObj)
        columns = signals["basic_signals"]

        print("start conversion")
        # Create a folder name "data" within the same dir of this script and place all the mat file in side the "data" folder.
        file_dict = {}

        for file in file_list:
            path = os.path.join(py_eval_path,"0_data_mat", file+".mat")
            file_dict[file] = path

        for file, path in file_dict.items():
            print(file)
            df = pd.DataFrame()
            converter = conversion(path, df)
            mat_file = converter.read_mat()
            df = converter.mat_to_df(mat_file, columns, remove_list)
            LC_raw = converter.laneChange_Cutter_raw()
            df = converter.laneChange_Maker(LC_raw)
            # All the labelled .csv file will be output and saved inside the "1_full_sequence_dataset" dir.
            converter.save_to_csv(file)

        cutter1 = cutter()
        cutter1.training_cutting(file_list)
        cutter1.sort_files('LaneChange')
        cutter1.sort_files('LeftRightTurn')
        cutter1.sort_files('Not')
        relabeler1 = relabeler()
        relabeler1.training_relabeling()
        os.rmdir(os.path.join(py_eval_path, '2_Cut_Data'))