import numpy as np
import pandas as pd
from operator import itemgetter
from itertools import groupby
import pathlib
import random
import os
import math
from utils.utilities import utilities
'''
This class exists to create the dataset appropriate for the Level of Dynamics analysis.
It will only work if the only csv files in the _temp folder are the converted measurement files.
To get these, make sure to comment out line 84 in utils/cleanup.py, save, run the evaluation of the measurement files that need to be analyzed and 
then copy the files from the "1_full_sequence_dataset" folder into the "_temp" folder.
'''
class dynamics_data_cutter():
    '''
    This method determines the ranges that are to be cut
    It always takes 1s before the blinker is set and 8s after the blinker is set
    '''
    def slicing(self, df, flag=False):
        data = df.index.tolist()
        ranges = []
        for k, g in groupby(enumerate(data),lambda x:x[0]-x[1]):
            group = (map(itemgetter(1), g))
            group = list(map(int, group))
            
            if flag:
                if group[-1]-100 > group[0]+21:
                    ranges.append((group[0], group[-1]))
            else:
                ranges.append((group[0], group[-1]))
        return ranges
    '''
    This method cuts the data based on the setting of the blinker
    If the blinker is set and the approve button is pushed it categorized as a LaneChange
    If the blinker is set and the approve button is not pushed it is categorized 
    as a LeftRightTurn (even if it is a roundabout or different corner case)
    '''
    def blinker_cutting(self, file_list, temp_path):
        util = utilities()
        py_eval_path = str(pathlib.Path(__file__).parent)
        drop_list = [
                        'GPS_time', 'q_T0', 'LaneChange_Approved',
                        'Left_Index', 'Right_Index', 'index',
                        'GPS_status',
                        'LatPos_abs', 'LongPos_abs',
                    ]
        cl = 0
        for file in file_list:
            if 'dynamic' in file:
                paths = [os.path.join(py_eval_path, 'LC_dataset_120', 'dynamic', file)]
                util.create_path(paths)
                outfolder = paths[0]
            elif 'comfort' in file:
                paths = [os.path.join(py_eval_path, 'LC_dataset_120', 'comfort', file),]
                util.create_path(paths)
                outfolder = paths[0]
            else:
                paths = [os.path.join(py_eval_path, 'LC_dataset_120', 'natural', file)]
                util.create_path(paths)
                outfolder = paths[0]
            df = pd.read_csv(os.path.join(temp_path, file+'.csv'))
            #the data is downsampled
            sample_rate = 10
            df = df.iloc[::sample_rate, :].copy()
            df = df.reset_index()
            #the ranges for the extracted segments are determined
            ranges = self.slicing(df.loc[(df['Left_Index'] == 1) | (df['Right_Index'] == 1)],  False)
            for index, _ in ranges:
                #LaneChange
                s1 = index-21
                e1 = index+100
                temp = util.remove_columns(df, drop_list)
                if temp[s1:e1].shape[0] > 120:
                    temp.iloc[s1:e1][:-1].to_csv(os.path.join(outfolder, str(cl)+'.csv'), index=False)
                    cl+=1

ddc = dynamics_data_cutter()
file_list = []
temp_path = os.path.join(str(pathlib.Path(__file__).resolve().parent.parent), '_temp')
for file in os.listdir(temp_path):
    if file[-4:] == '.csv':
        file_list.append(file[:-4])
ddc.blinker_cutting(file_list, temp_path)