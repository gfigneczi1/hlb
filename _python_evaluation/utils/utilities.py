import datetime
import pandas as pd
import numpy as np
import os
import math
from sklearn.preprocessing import MinMaxScaler
from scipy.signal import lfilter
from tensorflow.keras.utils import to_categorical
from operator import itemgetter
from itertools import groupby
import random
import datetime
'''
This class contains utility functions that are used at different points in the _python_evaluation folder
'''
class utilities():
    '''
    This method normalizes all of the data (between 0 and 1) for the classification model input
    '''
    def normalize(self, df):
        values = df.values.astype('float32')
        scaler = MinMaxScaler(feature_range=(0, 1))
        scaled = scaler.fit_transform(values)
        return scaled
    '''
    This method removes the any signals in the drop_list. Used to remove the signals not used for training
    '''
    def remove_columns(self, df, drop_list):
        temp = df.copy()
        for name in df:
            if name in drop_list:
                temp.drop(name, axis=1, inplace=True)
        temp = temp[temp.columns.drop(list(temp.filter(regex='Unnamed*')))]
        return temp
    '''
    This method filters all of the signals in the filter_list
    '''
    def data_filtering(self, df):
        filter_list = ['yawRateESP', 'AccelerationY_ESP', 'VelocityX_ESP', 'AccelerationX_ESP', ]
        n = 30  # the larger n is, the smoother curve will be
        b = [1.0 / n] * n
        a = 1
        filtered = df.copy()
        for column in df.columns:
            if column in filter_list:
                temp = df[column].to_numpy()
                filtered.loc[:, column] = lfilter(b, a, temp)
        return filtered
    '''
    This method filters a specific signal of choice
    '''
    def signal_filter(self, signal, n):
        # the larger n is, the smoother curve will be
        b = [1.0 / n] * n
        a = 1
        filtered_col = lfilter(b, a, signal)
        return filtered_col
    '''
    This method is used to create any needed path that doesn't exist yet
    '''
    def create_path(self, paths):
        for path in paths:
            if not os.path.exists(path):
                os.makedirs(path)
    '''
    This method calculates the relative acceleration accY_rel out of the yawRate the c2 value and the Velocity in x-direction
    '''
    def get_accY(self, df):
        c_2 = df['c2']
        v_x = df['VelocityX_ESP']
        yawR = df['yawRateESP']

        accY_rel = (yawR - c_2 * v_x) * v_x
        return accY_rel
    '''
    This function takes the left c01 and the right c01 of the data and returns the average c01
    '''
    def get_c01_average(self, df):
        c01_left = df['c01_left']
        c01_right = df['c01_right']
        c01_average = (c01_left + c01_right) / 2
        return c01_average
    '''
    This function finds the start and the end of the already labelled lane change column of choice
    '''
    def find_start_end(self, lc):
        prev = 0
        start = 0
        end = len(lc)
        for i in range(len(lc)):
            if prev == 0 and lc[i] == 1:
                start = i
                prev = 1
            if prev == 1 and lc[i] == 0:
                end = i
                break
        return start, end