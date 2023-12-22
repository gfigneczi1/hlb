import pandas as pd
import os
import numpy as np
import pathlib
from pathlib import Path
import tensorflow as tf
from tensorflow.keras.models import load_model
from numpy import split
from numpy import array
from operator import itemgetter
from itertools import groupby
from matplotlib import pyplot as plt
from _python_evaluation.utils.Relabel_LaneChange import relabeler
from _python_evaluation.utils.utilities import *
'''
This class is used to perform the inference seprate from training.
'''
class inference():
    def __init__(self, config):
        self.config = config
    
    '''
    This method determines the ranges that are to be cut
    It always takes 1s before the blinker is set and 8s after the blinker is set
    '''
    def slicing(self, df, flag=False):
        data = df.index.tolist()
        ranges = []
        for k, g in groupby(enumerate(data), lambda x:x[0]-x[1]):
            group = (map(itemgetter(1), g))
            group = list(map(int, group))
            if flag:
                if group[-1]-80 > group[0]+11:
                    ranges.append((group[0], group[-1]))
            else:
                ranges.append(group[0])
        return ranges
    '''
    This method loads the data for the classification. 
    The X_data consists of N files with 14 columns and 90 rows of data that is reshaped into shape (N, 36, 35, 1) 
    The y is simply a 1 if it is a lane change and a 0 if it is not
    It cuts the data and loads the data for the inference of the classification in one step.
    Two lists are synchronized to keep track of what start and end point one segment has within the longer sequence file. 
    '''
    def classification_data_loading(self, test_file, drop_list, sample_rate = 10, signals=12):
        X_data = []
        y = []
        #lists to keep track of the cut files and the (start,end) of the cut sequences in the long file
        LaneChange_files = []
        LaneChange_boundries = []
        cl = 0
        util = utilities()
    
        df = pd.DataFrame()
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        df = pd.read_csv(os.path.join(py_eval_path, '1_full_sequence_dataset', test_file+'.csv'))
        #downsampling
        df = df.iloc[::sample_rate, :].copy()
        df = df.reset_index()
        #determining ranges for slicing short sequences
        ranges = self.slicing(df.loc[(df['Left_Index'] == 1) | (df['Right_Index'] == 1)])
        #cutting similar to Extract_LaneChange file
        for index in ranges:
            s1 = index-11
            e1 = index+80
            temp = util.remove_columns(df, drop_list)
            if temp[s1:e1].shape[0] > 90:
                t = util.normalize(temp.iloc[s1:e1][:-1])
                if self.config == 'basic':
                    lc = t[:,-1]
                    temp.iloc[s1:e1][:-1].to_csv(os.path.join(py_eval_path, "2_Cut_Data", test_file, str(cl)+'.csv'), index=False)
                elif self.config == 'light':
                    temp.iloc[s1:e1][:-1].to_csv(os.path.join(py_eval_path, "2_Cut_Data", test_file, str(cl)+'.csv'), index=False)
                elif self.config == 'online':
                    lc = t[:,-2]
                    t = t[:,:-2]
                    temp.iloc[s1:e1][:-1].to_csv(os.path.join(py_eval_path, "2_Cut_Data", test_file, str(cl)+'.csv'), index=False)
                #X-data structured in the same way as during training
                if signals == 14:
                    input_shape = (36, 35, 1) 
                else:
                    input_shape = (36, 30, 1)
                X_data.append(t.reshape(input_shape))
                LaneChange_files.append(test_file + '/' + str(cl)+'.csv')
                LaneChange_boundries.append((s1, e1))
                cl += 1
                y.append(0)
        X_data = array(X_data)
        return X_data, y, LaneChange_files, LaneChange_boundries
    '''
    This method loads the data for the regression. 
    The X_data consists of N files with 2 columns (relative lateral acceleration, average c01) and 90 rows of data that is reshaped into shape (N, 90, 1, 2) 
    The y data is the final kinetic relabeling.
    '''
    def regression_data_loading(self, LaneChange_inferred):
        X_data = []
        y_data = []
        util = utilities()
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        for file in LaneChange_inferred:
            file_path = os.path.join(py_eval_path, '5_inferred', file)
            df = pd.read_csv(file_path)

            relabeler1 = relabeler()
            #calculating and loading X_data
            accY = util.get_accY(df)
            c01 = util.get_c01_average(df)
            out = np.column_stack((accY, c01))
            X_data.append(out)
            ##relabeling for y_data
            lc1, accY_start, accY_end = relabeler1.add_accY_labels(df)
            lc2, c01_start, c01_end = relabeler1.add_c01_labels(df)
            lc3 = relabeler1.kinetic_labels(lc1, accY_start, accY_end, lc2, c01_start, c01_end)
            start, end = util.find_start_end(lc3)
            df['kinetic_LaneChange'] = lc3
            df.to_csv(file_path, index=False)
            y_data.append((start, end))
        data_X = np.array(X_data)
        data_y = np.array(y_data)
        #reshape: (N, 90, 1, 2) with N samples, 90 in the time dimension, 1 data column per channel with 2 channels
        data_X = data_X.reshape(data_X.shape[0], 90, 1, 2)
        return data_X, data_y
    '''
    This method determines the direction of the lane change based on a c01 threshold for left and right each
    '''
    def get_direction(self, sub_path, util):
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        full_path = os.path.join(py_eval_path,'5_inferred', sub_path)
        df = pd.read_csv(full_path)
        c01_average = util.get_c01_average(df)
        direction = 0
        c01_thr_start = -0.5
        for i in range(len(c01_average)):
            if c01_average[i] < c01_thr_start:
                direction = 'left'
                break
            if c01_average[i] > c01_thr_start*(-1):
                direction = 'right'
                break
        return direction
    '''
    This method adds the prediction as an additional column to the data segments in the inferred folder
    '''
    def add_prediction(self, LaneChange_inferred, prediction):
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        for i, file in enumerate(LaneChange_inferred):
            file_path = os.path.join(py_eval_path, '5_inferred', file)
            df = pd.read_csv(file_path)
            pred_lc = []
            for j in range(90):
                pred_lc.append(0)
            pred_start = int(prediction[i, 0]) 
            pred_end = int(prediction[i, 1])
            if pred_start > 0 and pred_end < 91:
                for j in range(pred_start, pred_end):
                    pred_lc[j] = 1 
            df['prediction'] = pred_lc
            df.to_csv(file_path, index=False)
    '''
    This method performs the inference. The names for the models for classification and regression have to be input
    The 'test_file' input is a list of all of the files that are being inferred at once.
    The drop rate and sample rate are needed for the data loading and have to be input as well. 
    The 'signals' indicates the number of signals the classification model was trained with 
    '''
    def infer(self, class_model_name, regression_model_name, test_file, drop_list, sample_rate = 10, signals=12):
        util = utilities()
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        #creating paths to save the processed data in
        paths = [os.path.join(py_eval_path, '5_inferred', test_file), os.path.join(py_eval_path, '2_Cut_Data', test_file), os.path.join(py_eval_path,'6_results')]
        util.create_path(paths)
        '''CLASSIFICATION  '''
        #loading data for classification
        data, y, LaneChange_files, LaneChange_boundries = self.classification_data_loading(test_file, drop_list, sample_rate, signals=signals)
        if data.shape[0] == 0:
            print("No data to classify")
        else:
            #loading classification model
            class_model = load_model(py_eval_path + '/models/' + class_model_name + '.hdf5')
            #classification prediction
            yhat = class_model.predict(data,  verbose=0)
            #evaluation of classification prediction
            pre = []
            raw = []
            for i in range(len(yhat)):
                if yhat[i][0] > yhat[i][1]:
                    pre.append(1)
                else:
                    pre.append(0)
                raw.append(yhat[i])
            acc = 0
            temp = []
            acc = 0
            #list to keep track of all files that have been classified as Lane Changes
            #all the files that have not been classified as Lane Changes are discarded
            LaneChange_inferred = []
            LaneChange_boundries_inferred = []
            for x, y, z, a, b in zip(pre, raw, y, LaneChange_files, LaneChange_boundries):
                if x == 1:
                    os.replace(py_eval_path + '/2_Cut_Data/' + a, py_eval_path + '/5_inferred/' + a)
                    LaneChange_inferred.append(a)
                    LaneChange_boundries_inferred.append((b[0], b[1]))
                if x == z:
                    acc+=1
            print(yhat.shape[0], "files")
            LaneChange_boundries_inferred = array(LaneChange_boundries_inferred)
        if not os.listdir(py_eval_path +'/5_inferred/'+ test_file):
            print('no lane changes detected by classifier')
        else:
            '''REGRESSION'''
            #loading data for regression
            data_X, data_y = self.regression_data_loading(LaneChange_inferred)
            #loading model for regression
            regression_model = load_model(py_eval_path + '/models/' + regression_model_name + '.hdf5')
            #regression prediction
            prediction = regression_model.predict(data_X)
            #regression performance evaluation
            total_start = 0
            total_end = 0
            #list to keep track of the boundries found in the prediction of the classified lane changes
            #list will contain the start and end points of the lane changes in the long sequence and an indicator if it is a left or a right lane change
            Final_LaneChange_boundries = []
            for i in range(prediction.shape[0]):
                #finding direction for lane change
                direction = self.get_direction(LaneChange_inferred[i], util)
                #determining start and end in long sequence
                final_start = (LaneChange_boundries_inferred[i, 0]+prediction[i, 0])/10
                final_end = (LaneChange_boundries_inferred[i, 0]+prediction[i, 1])/10
                ground_truth_start = (LaneChange_boundries_inferred[i, 0]+data_y[i, 0])/10
                ground_truth_end = (LaneChange_boundries_inferred[i, 0]+data_y[i, 1])/10
                Final_LaneChange_boundries.append((ground_truth_start, ground_truth_end, final_start, final_end, direction))
                length_difference = (prediction[i,1]-prediction[i,0])-(data_y[i,1]-data_y[i,0])
                total_start += 10*abs(prediction[i,0] - data_y[i,0])
                total_end += 10*abs(prediction[i,1] - data_y[i,1])
            self.add_prediction(LaneChange_inferred, prediction)
            print('test mean absolute error: ', total_start/prediction.shape[0], total_end/prediction.shape[0])
            print('----------------------------')
            #if one file is inferred, the generated csv will include all of the starts and ends of the lane changes
            result = pd.DataFrame(array(Final_LaneChange_boundries))
            result.to_csv(py_eval_path + '/6_results/'+ test_file +'_results.csv', header=['kinetic_labels_starts', 'kinetic_labels_ends','regression_starts', 'regression_ends', 'direction']) #comment out in final implementation
            #saves to _temp folder
            top = py_eval_path = str(pathlib.Path(__file__).parent.parent.parent)
            save_path = os.path.join(top, '_temp', test_file+'_results.csv')
            result.to_csv(save_path, header=['kinetic_labels_starts', 'kinetic_label_ends','regression_starts', 'regression_ends', 'direction'])