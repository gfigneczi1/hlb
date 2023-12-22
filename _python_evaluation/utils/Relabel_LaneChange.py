import pandas as pd
import numpy as np
import os
import pathlib
from _python_evaluation.utils.utilities import *
'''
This class adds the kinetic labelling based on the relative lateral acceleration 
and the average of the left and right c01
The kinetic labelling is used as a ground truth for the training of the regression model
'''
class relabeler():
    '''
    This method creates the labelling based on the relative lateral acceleration
    It thresholds the relative lateral acceleration with 3 thresholds(start, mid, end) 
    and adds the settling time to the lane change to settle 
    '''
    def add_accY_labels(self, df):
        util = utilities()
        a_y_rel = util.get_accY(df) #relative lateral acceleration formula in utilities file
        a_y_rel = util.signal_filter(a_y_rel, 3)
        #three relative acceleration thresholds & settling time
        ACC_Y_trh_start = 0.009*35 
        ACC_Y_trh_mid = -0.0065*19
        ACC_Y_trh_end = -0.0054*36 
        settling_time = 14 
        
        lc1 = []
        for i in range(len(a_y_rel)):
            lc1.append(0)
        
        man_start_point = 0
        man_end_point = 0
        #loop input sequence
        for i in range(len(a_y_rel)):
            man_settled = False
            #first threshold check for left lane change
            if a_y_rel[i] > ACC_Y_trh_start:
                #start of maneuver set
                man_start_point = i
                for j in range(i, len(a_y_rel)):
                    #second theshold check left
                    if a_y_rel[j] < ACC_Y_trh_mid:
                        for k in range(j+10, len(a_y_rel)):
                            #third threshold check left
                            if a_y_rel[k] > ACC_Y_trh_end:
                                man_settled = True
                                #end of maneuver set
                                man_end_point = k
                                break
                    if man_settled:
                        #lane change time span is labelled with 1 
                        for s in range(man_start_point, man_end_point+settling_time):
                            if s+1 > len(a_y_rel):
                                break
                            lc1[s] = 1
                        break
            #first threshold check for right lane change
            if a_y_rel[i] < (ACC_Y_trh_start * (-1)):
                #start of maneuver set
                man_start_point = i
                for j in range(i, len(a_y_rel)):
                    #second threshold check right
                    if a_y_rel[j] > (ACC_Y_trh_mid*(-1)):
                        for k in range(j+10, len(a_y_rel)):
                            #third threshold check right
                            if a_y_rel[k] < (ACC_Y_trh_end*(-1)):
                                man_settled = True
                                #end of maneuver set
                                man_end_point = k
                                break
                    if man_settled:
                        #lane change time span is labelled with 1
                        for s in range(man_start_point, man_end_point+settling_time):
                            if s+1 > len(a_y_rel):
                                break
                            lc1[s] = 1
                        break
            if man_settled:
                break
        #to fix a bug with below formula that would wrongly return a label_end time of 14 for non lane changes
        if man_end_point == 0:
            settling_time = 0
        #final end of the lane change set
        label_end = man_end_point+settling_time
        return lc1, man_start_point, label_end
    '''
    This method creates the labelling based on the relative lateral acceleration
    It thresholds the average of c01_left and c01_right with 4 thresholds(start, mid1, mid2, end) 
    and adds the settling time to the lane change to settle 
    '''
    def add_c01_labels(self, df):
        util = utilities()
        #c01 formula in utility function
        c01_average = util.get_c01_average(df)
        #4 thresholds for average c01
        c01_thr_start = -0.5
        c01_thr_mid1 = -1.2
        c01_thr_mid2 = 1 #change to 1
        c01_thr_end = 0.5
        settling_time = 14
        
        lc2 = []
        for i in range(len(c01_average)):
            lc2.append(0)
        man_start_point = 0
        man_end_point = 0
        #loop input sequence
        for i in range(len(c01_average)):
            man_settled = False
            #first threshold check for left lane change
            if c01_average[i] < c01_thr_start:
                #start of maneuver set
                man_start_point = i
                for j in range(i, len(c01_average)):
                    #second threshold check left
                    if c01_average[j] < c01_thr_mid1:
                        for l in range(j, len(c01_average)):
                            #third threshold check left
                            if c01_average[l] > c01_thr_mid2:
                                for k in range(l+5, len(c01_average)):
                                    #fourth threshold check left
                                    if c01_average[k] < c01_thr_end:
                                        man_settled = True
                                        #end of maneuver set
                                        man_end_point = k
                                        break
                    if man_settled:
                        #lane change time span is labelled with 1
                        for s in range(man_start_point, man_end_point+settling_time):
                            if s+1 > len(c01_average):
                                break
                            lc2[s] = 1 
                        break
            #first threshold check for right lane change
            if c01_average[i] > c01_thr_start*(-1):
                #start of maneuver set
                man_start_point = i
                for j in range(i, len(c01_average)):
                    #second threshold check right
                    if c01_average[j] > c01_thr_mid1*(-1):
                        for l in range(j, len(c01_average)):
                            #third threshold check right
                            if c01_average[l] < c01_thr_mid2*(-1):
                                for k in range(l+5, len(c01_average)):
                                    #fourth threshold check right
                                    if c01_average[k] > c01_thr_end*(-1):
                                        man_settled = True
                                        #end of maneuver set
                                        man_end_point = k
                                        break
                            if man_settled:
                                break
                    if man_settled:
                        #lane change time span is labelled with 1
                        for s in range(man_start_point, man_end_point+settling_time):
                            if s+1 > len(c01_average):
                                break
                            lc2[s] = 1
                        break
            if man_settled:
                break
        #to fix a bug with below formula that would wrongly return a label_end time of 14 for non lane changes
        if man_end_point == 0:
            settling_time = 0
        #final end of the lane change set
        label_end = man_end_point+settling_time
        return lc2, man_start_point, label_end
    '''
    This method creates the final labelling used as the ground truth for the training of the regression model
    In general the start of the relative lateral acceleration label and the end of the average c01 label are chosen 
    as the bounds of the lane change
    '''        
    def kinetic_labels(self, lc1, accY_start, accY_end, lc2, c01_start, c01_end):
        lc3 = []
        for i in range(len(lc1)):
            lc3.append(0)
        #check if c01_end is within the bounds of the file  
        if c01_end > len(lc2):
            c01_end = len(lc2)
        #check if the lateral acceleration label start is within a certain time range before the c01 label starts
        #c01 label is in general delayed compared to the acceleration label -> acceleration label returns the preferred start
        #if the difference of starts between acceleration and c01 label is too large, then the acceleration is triggered by a disturbance
        if accY_start in range(c01_start-20, c01_start):
            for i in range(accY_start, c01_end):
                lc3[i] = 1
        #if the lane change is still detected by the c01 label despite the difference of the starts 
        #we take 1s ahead of the c01 label start as the final kinetic label start
        if c01_end > 0:
            if c01_start > 10:
                for i in range(c01_start-10, c01_end):
                    lc3[i] = 1
            else:
                for i in range(c01_start, c01_end):
                    lc3[i] = 1
        return lc3
    '''
    This method does the relabeling during the training process
    It is called in the preprocessing script
    '''            
    def training_relabeling(self):
        print("start relabeling")
        util = utilities()
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        #iterate through categories in 3_cut_dataset: train, val, test
        for cat in os.listdir(os.path.join(py_eval_path, "3_cut_dataset")):
            print(cat)
            #iterate through folders: LaneChange, LeftRightTurn, Not
            for folder in os.listdir(os.path.join(py_eval_path, "3_cut_dataset", cat)):
                #iterate through individually cut up lane change files
                for file in os.listdir(os.path.join(py_eval_path, "3_cut_dataset", cat, folder)):
                    dir = os.path.join(py_eval_path, "3_cut_dataset", cat, folder, file)
                    df = pd.read_csv(dir)
                    #add the relative lateral acceleration labeling
                    lc1, accY_start, accY_end = self.add_accY_labels(df)
                    df['accY_LaneChange'] = lc1
                    #add the average c01 labeling
                    lc2, c01_start, c01_end = self.add_c01_labels(df)
                    df["c01_LaneChange"] = lc2
                    #add the final labelling
                    lc3 = self.kinetic_labels(lc1, accY_start, accY_end, lc2, c01_start, c01_end)
                    df['kinetic_LaneChange'] = lc3
                    #add relabeled file to realabeled_dataset
                    out_path = [os.path.join(py_eval_path, "4_relabeled_dataset", cat , folder)]
                    util.create_path(out_path)
                    df.to_csv(os.path.join(py_eval_path, '4_relabeled_dataset', cat, folder, file), index=False)