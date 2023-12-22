import pandas as pd
import numpy as np
from dynamics_utils.data_getter import data_getter
from utils.utilities import utilities

class dynamics_segmentation():
    def segment_lane_change(self, a_y_rel, df):
        thr1, thr2, thr3, thr4 = self.segment_c01_thr(df)
        seg1 = a_y_rel[0:thr1]
        seg2 = a_y_rel[thr1:thr4]
        direction = self.get_direction(a_y_rel)
        if len(seg1) == 0 or len(seg2) == 0:
            t1 = thr1
            t2 = thr2
        elif direction == 'left':
            t1 = self.get_max_peak_time(seg1)+5
            t1 = t1
            t2 = self.get_min_peak_time(seg2)+5
            t2 = thr1 + t2
        elif direction == 'right':
            t1 = self.get_min_peak_time(seg1)+5
            t1 = t1
            t2 = self.get_max_peak_time(seg2)+5
            t2 = thr1 + t2
        else:
            t1 = 0
            t2 = 0 
        return t1, t2
    def segment_alternate(self, a_y_rel, df):
        t1, t2 = 0, 0 
        midpoint = 3
        util = utilities()
        c01 = util.get_c01_average(df)
        if len(c01) == 0 or len(a_y_rel) == 0:
            return 30, 60
        max_c01 = max(c01)
        for i, element in enumerate(c01):
            if element == max_c01:
                midpoint = i
        if midpoint < 3:
            midpoint = 3
        if c01[midpoint-3] < 0:
            direction = 'left'
        else:
            direction = 'right'
        seg1 = a_y_rel[0:midpoint]
        seg2 = a_y_rel[midpoint:len(a_y_rel)]
        if len(seg1) == 0 or len(seg2) == 0:
            return 0, 0
        elif direction == 'left':
            t1 = self.get_max_peak_time(seg1)+2
            t2 = self.get_min_peak_time(seg2)+2
            t2 = midpoint + t2
        elif direction == 'right':
            t1 = self.get_min_peak_time(seg1)+2
            t2 = self.get_max_peak_time(seg2)+2
            t2 = midpoint + t2
        else:
            t1 = 0
            t2 = 0
        return t1, t2
    def get_max_peak_time(self, a_y_rel):
        t = 0
        peak = max(a_y_rel)
        for i, a in enumerate(a_y_rel):
            if peak == a:
                t = i
                break
        return t
    def get_min_peak_time(self, a_y_rel):
        t = 0
        peak = min(a_y_rel)
        for i, a in enumerate(a_y_rel):
            if peak == a:
                t = i
                break
        return t
        
    def get_direction(self, a_y_rel):
        ACC_Y_trh_start = 0.009*35 
        for i in range(len(a_y_rel)):
            if a_y_rel[i] < ACC_Y_trh_start:
                direction = 'left'
                break
            if a_y_rel[i] > ACC_Y_trh_start*(-1):
                direction = 'right'
                break
        return direction
    '''LANE CHANGE SEGMENTATION'''
    def segment_ay_thr(self, a_y_rel):
        t1 = 0
        t2 = 0
        t3 = 0
        t4 = 0
        #thresholds
        ACC_Y_trh_start = 0.009*35 
        ACC_Y_trh_mid = -0.0065*19
        ACC_Y_trh_end = -0.0054*36 
        settling_time = 20#14
        #loops to check thresholds
        for i in range(len(a_y_rel)):
            if a_y_rel[i] > ACC_Y_trh_start:
                t2 = i
                for j in range(i, len(a_y_rel)):
                    if a_y_rel[j] < ACC_Y_trh_mid:
                        t3 = j
                        for k in range(j+10, len(a_y_rel)):
                            #third threshold check left
                            if a_y_rel[k] > ACC_Y_trh_end:
                                t4 = k
                                break
                    if t4 > 0:
                        break
            if a_y_rel[i] < (ACC_Y_trh_start * (-1)):
                t2 = i
                for j in range(i, len(a_y_rel)):
                    if a_y_rel[j] > (ACC_Y_trh_mid*(-1)):
                        t3=j
                        for k in range(j+10, len(a_y_rel)):
                            #third threshold check right
                            if a_y_rel[k] < (ACC_Y_trh_end*(-1)):
                                t4 = k
                                break
                    if t4 > 0:
                        break
            if t2 > 0:
                break
        if t4 > 0:
            if t4 + settling_time < len(a_y_rel):
                t4 = t4 + settling_time
            else:
                t4 = len(a_y_rel)
        return t1, t2, t3, t4

    def segment_c01_thr(self, df):
        t1, t2, t3, t4 = 0, 0, 0, 0
        util = utilities()
        #c01 formula in utility function
        c01_average = util.get_c01_average(df)
        #4 thresholds for average c01
        c01_thr_start = -0.5
        c01_thr_mid1 = -1.2
        c01_thr_mid2 = 1 #change to 1
        c01_thr_end = 0.5
        settling_time = 14
        #loop input sequence
        for i in range(len(c01_average)):
            man_settled = False
            #first threshold check for left lane change
            if c01_average[i] < c01_thr_start:
                t1 = i
                for j in range(i, len(c01_average)):
                    #second threshold check left
                    if c01_average[j] < c01_thr_mid1:
                        t2 = j
                        for l in range(j, len(c01_average)):
                            #third threshold check left
                            if c01_average[l] > c01_thr_mid2:
                                t3 = l
                                for k in range(l+5, len(c01_average)):
                                    #fourth threshold check left
                                    if c01_average[k] < c01_thr_end:
                                        man_settled = True
                                        t4 = k
                                        break
                    if man_settled:
                        break
            #first threshold check for right lane change
            if c01_average[i] > c01_thr_start*(-1):
                t1 = i
                for j in range(i, len(c01_average)):
                    #second threshold check right
                    if c01_average[j] > c01_thr_mid1*(-1):
                        t2 = j
                        for l in range(j, len(c01_average)):
                            #third threshold check right
                            if c01_average[l] < c01_thr_mid2*(-1):
                                t3 = l
                                for k in range(l+5, len(c01_average)):
                                    #fourth threshold check right
                                    if c01_average[k] > c01_thr_end*(-1):
                                        man_settled = True
                                        t4 = k
                                        break
                            if man_settled:
                                break
                    if man_settled:
                        break
            if man_settled:
                break
        #to fix a bug with below formula that would wrongly return a label_end time of 14 for non lane changes
        if t4 == 0:
            settling_time = 0
        #final end of the lane change set
        t4 = t4+settling_time
        return t1, t2, t3, t4