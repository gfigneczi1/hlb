import pandas as pd
import numpy as np
import datetime
import os
import pathlib
from _python_evaluation.utils.utilities import utilities

class automatic_evaluation():
    
    def __init__(self, gentle_time_threshold=100, strict_time_threshold=30):
        self.gentle_time_threshold = gentle_time_threshold
        self.strict_time_threshold = strict_time_threshold
    '''
    This method checks if one the labels cannot detect any lane change. It is used to compare the kinetic label to the regression prediction.
    If the ground truth and the prediction detect no lane change it is a True Negative.
    If the ground truth fails and the prediction has a lane change it is a False Positive.
    If the ground truth detects a lane change and the prediction fails to do so, it is a False Negative.
    Evaluates one row at a time
    '''
    def basic_check(self, TN, FP, FN, a, b, c, d, causes):
        #the ground truth is negative
        if b == a+900:
            #ground truth and prediction are negative
            if d >= c+899 and d <=c+901:
                TN += 1
            #the prediction is positive; only ground truth is negative
            else:
                FP += 1
                causes[0, 0] += 1
        #only prediction is negative
        elif d >= c+899 and d <=c+901:
            FN += 1
            causes[1,0] += 1
        return TN, FP, FN, causes
    '''
    Checks by how much the prediction start and end times deviate from the ground truth start and end times.
    Two thresholds are used to compare and then categorize.
    Since the early start and the late end can deviate more than the late start and the early end, a 'gentle threshold' is used to allow for more leeway compared to the strict one. 
    Evaluates one row at a time.
    '''
    def time_compare(self, y_start, y_end, yhat_start, yhat_end):
        sd = y_start - yhat_start #start distance
        ed = y_end - yhat_end #end distance
        start_eval = 'x'
        end_eval = 'x'
        if sd <= self.strict_time_threshold and sd >= -self.strict_time_threshold:
            start_eval = 'sg' #sg = start good
        elif sd > self.strict_time_threshold and sd <= self.gentle_time_threshold:
            start_eval = 'se1' #se = start early fine
        elif sd > self.strict_time_threshold and sd > self.gentle_time_threshold:
            start_eval = 'se2' #se = start early bad
        elif sd < -self.strict_time_threshold:
            start_eval = 'sl' #sl = start late
        if ed <= self.strict_time_threshold and ed >= -self.strict_time_threshold:
            end_eval = 'eg' #eg = end good
        elif ed > self.strict_time_threshold:
            end_eval = 'ee' #ee = end early
        elif ed < -self.strict_time_threshold and ed >= -self.gentle_time_threshold:
            end_eval = 'el1' #el = end late fine
        elif ed < -self.strict_time_threshold and ed < -self.gentle_time_threshold:
            end_eval = 'el2' #el = end late bad
        return start_eval, end_eval
    '''
    Based on the category the lane change is evaluated:
    It is a true positive if the start is either good or slightly early and the end is good or slightly late
    It is a false positive if the start is far too early or the end is far too late.
    It is a false negative if the start is too late or the end is too early.
    Evaluates one row at a time.
    '''
    def time_evaluate(self, start_eval, end_eval, TP, FP, FN, causes):
        if (start_eval == 'sg' or start_eval == 'se1') and (end_eval == 'eg' or end_eval == 'el1'):
            TP += 1
        elif (start_eval == 'se2' and end_eval == 'eg'):
            FP += 1
            causes[0,1] += 1
        elif (start_eval == 'sg' and end_eval == 'el2'):
            FP += 1
            causes[0,2] += 1
        elif (start_eval == 'se2' and end_eval == 'el2'):
            FP += 1
            causes[0,1] += 1
            causes[0,2] += 1
        elif ((start_eval == 'se1' or start_eval == 'se2') and end_eval == 'ee') or (start_eval == 'sl' and (end_eval == 'el1' or end_eval == 'el2')):
            FN += 1
            causes[1,1] += 1
            causes[1,2] += 1
        elif (start_eval == 'sl' and end_eval == 'ee'):
            FN += 1
            causes[1,1] += 1
            causes[1,2] += 1
        elif (start_eval == 'sg' and end_eval == 'ee'):
            FN += 1
            causes[1,2] += 1
        elif (start_eval == 'sl' and end_eval == 'eg'):
            FN += 1
            causes[1,1] += 1
        return TP, FP, FN, causes
    '''
    Currently the ground truth data is the kinetic (relative lateral acceleration & average c01) labels and the predictions are those of the regression model 
    When this method is called, the evaluation matrix comparing the results of the regression model and the results of the kinetic label as the ground truth is produced.
    It saves the evaluation matrix to _temp.
    '''
    def evaluate(self, ground_truth='kinetic_labels', prediction='regression'):
        #counters
        n = 1
        TP = 0
        FP = 0
        FN = 0
        TN = 0
        causes = np.array([[0, 0, 0, 0],[0, 0, 0, 0]]) #[[FP_basic, FP_start, FP_end, FP_both],[FN_basic, FN_start, FN_end, FN_both]]
        #import results
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        util = utilities()
        for file in os.listdir(py_eval_path + '/6_results'):
            if not 'online' in file: #should not be needed if 6_results is cleaned after every evaluation call
                df = pd.read_csv(os.path.join(py_eval_path, '6_results',file))
                y_starts = df[ground_truth + '_starts']*100
                y_ends = df[ground_truth + '_ends']*100
                yhat_starts = df[prediction + '_starts']*100 
                yhat_ends = df[prediction + '_ends']*100
                for a, b, c, d in zip(y_starts, y_ends, yhat_starts, yhat_ends):
                    n += 1
                    if b == a + 900 or (d >= c+899 and d <=c+901):
                        TN, FP, FN, causes = self.basic_check(TN, FP, FN, a, b, c, d, causes)
                    else:
                        start_eval, end_eval = self.time_compare(a, b, c, d)
                        TP, FP, FN, causes = self.time_evaluate(start_eval, end_eval, TP, FP, FN, causes)
        evaluation_matrix = [['TP: ' + str(TP), 'FP: ' + str(FP)],['FN: ' + str(FN), 'TN: '+ str(TN)]]
        evaluation_matrix = pd.DataFrame(evaluation_matrix)
        print('number of files:', n-1)
        print('evaluation matrix')
        print(evaluation_matrix)
        print("Failed to label(positive pred.):",causes[0, 0])
        print("Failed to label(negative pred.):",causes[1, 0])
        causes[0, 3] = (causes[0, 0]+causes[0, 1]+causes[0, 2])-FP
        causes[1,3] = (causes[1, 0]+causes[1, 1]+causes[1,2])-FN
        print("FP caused by start:  ", causes[0, 1]-causes[0,3])
        print("FN caused by start:  ",causes[1,1]-causes[1,3])
        print("FP caused by end:    ",causes[0,2]-causes[0,3])
        print("FN caused by end:    ",causes[1,2]-causes[1,3])
        print("FP caused by both:    ",causes[0,3])
        print("FN caused by both:    ",causes[1,3])
        if n > 1:
            print('Accuracy: ', (TP+TN)/(n-1)) #how many of the files were predicted correctly
        top = str(pathlib.Path(__file__).resolve().parent.parent.parent)
        evaluation_matrix.to_csv(os.path.join(top, '_temp', 'Evaluation_matrix_'+ground_truth+'_'+prediction+'.csv'))
    '''
    Makes the evaluation matrix with the ground truth (can be the kinetic labels or the regression model) and the online prediction.
    It checks whether the start and the end are in close enough to the ground truth, which can be either the results of the kinetic labeling or the results of the regression.
    pred_start has to be in [ground_truth_start-gentle_thr, ground_truth_start+strict_thr]
    pred_end has to be in [ground_truth_start-strict_thr, ground_truth_start+gentle_thr]
    For every True Positive, the online prediction of the direction of the lane change is compared with the ground truths direction prediction of the lane change.
    '''
    def compare_online(self, ground_truth_starts, ground_truth_ends, starts, ends, ground_truth_directions, directions):
        FP = 0
        TP = 0
        FN = 0
        n = len(ground_truth_starts)
        n_lc = np.arange(0, n)
        m = np.arange(0, len(starts))
        True_Positives = []
        for a, b, j in zip(starts, ends, m):
            #reads out the values that can be found in 6_results files
            for c, d, f, i in zip(ground_truth_starts, ground_truth_ends, ground_truth_directions, n_lc):
                if a in range(int(c-self.gentle_time_threshold),int(c+self.strict_time_threshold)) and b in range(int(d-self.strict_time_threshold), int(d+self.gentle_time_threshold)):
                    TP += 1
                    #The online method predicts a left lane change
                    if directions[b] == 1:
                        True_Positives.append([f, int(c/100), int(d/100), 'left', a/100, b/100, f == 'left'])
                    #The online method predicts a right lane change
                    elif directions[b] == 2:
                        True_Positives.append([f, int(c/100), int(d/100), 'right', a/100, b/100, f == 'right'])
                    #It is unknown what direction the lane change takes (0 in online method)
                    else:
                        True_Positives.append([f, int(c/100), int(d/100), directions[j], a/100, b/100, False])
                    break
                elif i == n-1:
                    FP += 1
        if TP < n:
            FN = n - TP
        eval_matrix = [['TP='+str(TP), 'FP='+str(FP)],['FN='+str(FN), 'TN=0']]
        eval_matrix = pd.DataFrame(eval_matrix)
        True_Positives = pd.DataFrame(np.array(True_Positives))
        return eval_matrix, True_Positives
    '''
    This method evaluates the online lane change prediction and is called during the python evaluation. It saves the evaluation matrix to _temp.
    '''
    def evaluate_online(self, file, starts, ends, directions):
        print('online evaluation of:', file)

        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        top = str(pathlib.Path(__file__).parent.parent.parent)
        util = utilities()
        
        file_path = py_eval_path + '/6_results/' + file + '_results.csv'
        df = pd.read_csv(file_path)
        kinetic_labels_starts = df['kinetic_labels_starts']*100
        kinetic_labels_ends = df['kinetic_labels_ends']*100
        ground_truth_directions = df['direction']
        eval_matrix_kinetic_gt, kinetic_True_Positives = self.compare_online(kinetic_labels_starts, kinetic_labels_ends, starts, ends, ground_truth_directions, directions)
        
        regression_starts = df['regression_starts']*100
        regression_ends = df['regression_ends']*100
        eval_matrix_reg_gt, regression_True_Positives = self.compare_online(regression_starts, regression_ends, starts, ends, ground_truth_directions, directions)
        
        print('evaluation matrix with the kinetic label as ground truth:')
        print(eval_matrix_kinetic_gt)
        eval_matrix_kinetic_gt.to_csv(os.path.join(top, '_temp', file+'_evaluation_matrix_kinetic_online.csv'))
        print('evaluation matrix with the regression prediction as ground truth:')
        print(eval_matrix_reg_gt)
        eval_matrix_reg_gt.to_csv(os.path.join(top, '_temp', file+'_evaluation_matrix_regression_online.csv'))
        
        if kinetic_True_Positives.shape[0] > 0:
            print("True positives with label as ground truth: ")
            print(kinetic_True_Positives)
            kinetic_True_Positives.to_csv(os.path.join(top, '_temp', file+'_True_Positives_Kinetic_labels_gt.csv'), header=['direction_gt','start_gt','end_gt','direction', 'start', 'end', 'direction_correct?'])
        if regression_True_Positives.shape[0] > 0:
            print("True positives with regression as ground truth: ")
            print(regression_True_Positives)
            regression_True_Positives.to_csv(os.path.join(top, '_temp', file+'_True_Positives_regression_gt.csv'), header=['direction_gt','start_gt','end_gt','direction', 'start', 'end', 'direction_correct?'])