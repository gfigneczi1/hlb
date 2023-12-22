import pathlib
import os
import numpy as np
import pandas as pd
import scipy
from scipy.signal import lfilter
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from utils.utilities import utilities
class data_getter():

    def derive(self, f, h):
        dx = []
        for i in range(len(f)-h):
            dx.append((f[i+h]-f[i])/h)
        for i in range(h):
            dx.append(dx[len(f)-h-1])
        return dx
    def column_filter(self, col, n):
        # the larger n is, the smoother curve will be
        b = [1.0 / n] * n
        a = 1
        filtered_col = lfilter(b, a, col)
        return filtered_col
    '''TO GET DATA'''
    def get_continuous_data(self, path):
        df = pd.read_csv(path)
        util = utilities()
        #acceleration y
        a_y_rel = util.get_accY(df)
        a_y_rel = self.column_filter(a_y_rel, 5)
        jerk = self.derive(a_y_rel, 1)
        jerk = np.array(self.column_filter(jerk, 5))
        return a_y_rel, jerk
    '''INIT JERK in range [time1, time2]'''
    def get_init_jerk(self, jerk, time1, time2):
        pos_peaks, _ = find_peaks(jerk[time1:time2])
        neg_peaks, _ = find_peaks((-1)*jerk[time1:time2])
        jerk_peaks = []
        for peak in pos_peaks:
            jerk_peaks.append(jerk[peak])
        for peak in neg_peaks:
            jerk_peaks.append(jerk[peak])
        if not jerk_peaks:
            return 0
        else:
            init_jerk = max(jerk_peaks)
            return init_jerk
    '''Maximum acceleration in range [time1, time2]'''
    def get_max_a_y(self, a_y_rel, time1, time2):
        peaks, _ = find_peaks(a_y_rel[time1:time2], height=0)
        min_peaks, _ = find_peaks(((-1)*a_y_rel[time1:time2]), height=0)
        a_y_peaks = []
        for peak in peaks:
            a_y_peaks.append(a_y_rel[peak])
        for peak in min_peaks:
            a_y_peaks.append((-1)*a_y_rel[peak])
        if not a_y_peaks:
            return 0
        else:
            max_a_y = max(a_y_peaks)
            return max_a_y
    '''GET INTEGRAL JERK of range [time1, time2]'''
    def get_integral_jerk(self, jerk, time1, time2):
        time_diff = time2 - time1
        if time_diff > 1:
            sq_jerk = np.transpose(jerk)*jerk
            sq_jerk_sub = sq_jerk[time1:time2]
            sq_jerk_sub = np.array(self.column_filter(sq_jerk_sub, 1))
            jerk_int = (1/time_diff)*np.sum(sq_jerk_sub)
            return jerk_int
        else:
            return 0
