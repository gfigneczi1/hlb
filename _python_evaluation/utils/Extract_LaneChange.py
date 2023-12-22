import numpy as np
import pandas as pd
from operator import itemgetter
from itertools import groupby
import pathlib
import random
import os
import math
from _python_evaluation.utils.utilities import *
'''
This class contains methods that cut the longer sequences from the 1_full_sequence_dataset into short files with 9s worth of data for training
and save them to the 2_Cut_Data folder. Another function then saves the reordered files to 3_cut_dataset
'''
class cutter():
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
                if group[-1]-80 > group[0]+11:
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
    def blinker_cutting(self, file_list, drop_list, paths):
        util = utilities()
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        cl = 0
        clr = 0
        for file in file_list:
            df = pd.read_csv(os.path.join(py_eval_path, '1_full_sequence_dataset', file+'.csv'))
            #the data is downsampled
            sample_rate = 10
            df = df.iloc[::sample_rate, :].copy()
            df = df.reset_index()
            #the ranges for the extracted segments are determined
            ranges = self.slicing(df.loc[(df['Left_Index'] == 1) | (df['Right_Index'] == 1)],  False)
            for index, _ in ranges:
                #LaneChange
                if df['LaneChange'][index] == 1:
                    s1 = index-11
                    e1 = index+80
                    temp = util.remove_columns(df, drop_list)
                    if temp[s1:e1].shape[0] > 90:
                        temp.iloc[s1:e1][:-1].to_csv(os.path.join(paths[0], str(cl)+'.csv'), index=False)
                        cl+=1
                #LeftRightTurn
                if df['LaneChange'][index] == 0:
                    s2 = index-11
                    e2 = index+80
                    temp = util.remove_columns(df, drop_list)
                    if temp[s2:e2].shape[0] > 90:
                        temp.iloc[s2:e2][:-1].to_csv(os.path.join(paths[1], str(clr)+'.csv'), index=False)
                        clr+=1
    '''
    This method cuts random parts of the long sequences where no blinker is set
    It basically generates non-LaneChange files for the classification model to not suffer
    from the drawbacks of data imbalance 
    '''
    def non_lane_change_cutting(self, file_list, drop_list, paths):
        util = utilities()
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        t = 0
        for file in file_list:
            tdf = pd.read_csv(py_eval_path+'/1_full_sequence_dataset/'+file+'.csv')
            sample_rate = 10
            tdf = tdf.iloc[::sample_rate, :].copy()
            tdf = tdf.reset_index()
            ranges = self.slicing(tdf.loc[tdf['LaneChange'] != 1], True)
            breaker = False
            for s, e in ranges:
                for i in range(2):
                    nbr = random.randint(s+11, e-80) 
                    s1 = nbr-11
                    e1 = nbr+80
                    temp = util.remove_columns(tdf, drop_list)
                    if temp[s1:e1].shape[0] > 90:
                        temp.iloc[s1:e1][:-1].to_csv(os.path.join(paths[2], str(t)+'.csv'), index=False)
                        t+=1
                    
                    if t == 400:
                        breaker = True
                        break
                if breaker:
                    break
            if breaker:
                break
    '''
    The method is executed in preprocessing to initiate the cutting process
    '''    
    def training_cutting(self, file_list):
        print("start cutting")
        #This is a list of the data signals that should not be used for training
        drop_list = [
                        'GPS_time', 'q_T0', 'LaneChange_Approved',
                        'Left_Index', 'Right_Index', 'index',
                        'GPS_status',
                        'LatPos_abs', 'LongPos_abs',
                    ]
        util = utilities()
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        paths = [os.path.join(py_eval_path, '2_Cut_Data','LaneChange'), os.path.join(py_eval_path,'2_Cut_Data','LeftRightTurn'), os.path.join(py_eval_path,'2_Cut_Data','Not')]  
        util.create_path(paths)
        self.blinker_cutting(file_list, drop_list, paths)
        self.non_lane_change_cutting(file_list, drop_list, paths)
    '''
    This method exists to sort the cut files from 2_Cut_Data into train, val and test in 3_cut_dataset
    The distribution is: 85% train, 10% validation, 5% test
    It is also executed during the preprocessing
    '''
    def sort_files(self, type):
        util = utilities()
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        in_path = os.path.join(py_eval_path, '2_Cut_Data', type)
        file_list = []
        for file in os.listdir(in_path):
            file_list.append(file)
        random.shuffle(file_list)
        paths = [os.path.join(py_eval_path, '3_cut_dataset', 'train', type), os.path.join(py_eval_path,'3_cut_dataset', 'val', type), os.path.join(py_eval_path, '3_cut_dataset', 'test', type)]
        util.create_path(paths)
        #Set split here: 85% train, 10% val, 5% test
        #new split: 70%, 15%, 15%
        train_end = math.ceil(0.7 * len(file_list))
        val_end = train_end + math.ceil(0.15 * len(file_list))
        for i in range(len(file_list) ):
            old_path = os.path.join(in_path, file_list[i])
            if i <= train_end:
                new_path = os.path.join(paths[0], file_list[i])
            elif i <= val_end and i > train_end:
                new_path = os.path.join(paths[1], file_list[i])
            elif i > val_end:
                new_path = os.path.join(paths[2], file_list[i])
            os.replace(old_path, new_path)
        os.rmdir(in_path)