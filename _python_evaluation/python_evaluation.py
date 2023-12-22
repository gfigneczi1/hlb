import os
import pandas as pd
import pathlib
import numpy
import math
import json
from tkinter import * 
import tkinter.messagebox 
from _python_evaluation.utils.Convert_LaneChange import conversion
from _python_evaluation.utils.Infer_LaneChange import inference
from _python_evaluation.utils.automatic_evaluation import automatic_evaluation
from _python_evaluation.utils.file_getter import file_getter
from _python_evaluation.utils.utilities import utilities
from _python_evaluation.utils.cleanup import cleanup

class python_evaluation():
    '''The class is initialized with the Python configuration (Python Evaluation profile) that is used to evaluate and the two threshold used for the automatic evaluation.'''
    def __init__(self, config, strict_time_threshold, gentle_time_threshold):
        self.config = config
        self.strict_time_threshold = int(100*float(strict_time_threshold))
        self.gentle_time_threshold = int(100*float(gentle_time_threshold))
    '''
    This method is for the online prediction evaluation and finds the starts and the ends of the online predictions.
    The evaluated online evaluation marks lane changes differently. If it thinks that a lane change may have started the elements in the columns below that point are 1.
    If it then confirmes that the lane change did indeed happen, the a few elements below that are 2.
    If it goes from 1 to 0 directly that means that the initially marked lane change was not a lane change.
    '''
    def find_online_starts_ends(self, lc):
        prev = 0
        potential_starts = []
        starts = []
        ends = []
        for i, element in enumerate(lc):
            if element == 0:
                prev = 0
            if element == 1 and prev == 0:
                prev = 1
                potential_starts.append(i)
            elif element == 1:
                prev = 1
            if element == 2 and prev == 1:
                prev = 2
                start = potential_starts.pop()
                starts.append(start)
                ends.append(i)
            elif element == 2:
                prev = 2
        print(starts)
        print(ends)
        return starts, ends
    '''
    This method executes the python_evaluation.
    '''
    def execute(self, class_model, reg_model):
        print('Python evaluation started')
        util = utilities()
        py_eval_path = str(pathlib.Path(__file__).parent)
        top = str(pathlib.Path(__file__).parent.parent)
        class_model_path = os.path.join(str(pathlib.Path(__file__).parent), 'models', class_model+'.hdf5')
        reg_model_path = os.path.join(str(pathlib.Path(__file__).parent), 'models', reg_model+'.hdf5')
        #create paths needed for python evaluation
        paths = [os.path.join(py_eval_path, '0_data_mat'), os.path.join(py_eval_path,'1_full_sequence_dataset'), os.path.join(py_eval_path,'6_results')]
        util.create_path(paths)
        #check if the classification model exists and stop execution with warning popup if it does not
        if not os.path.isfile(class_model_path):
            root=Tk() 
            root.withdraw()
            info = tkinter.messagebox.showinfo("ATTENTION","This Classification model does not exist!!")
            root.destroy()
        #check if the classification model exists and stop execution with warning popup if it does not
        elif not os.path.isfile(reg_model_path):
            root=Tk() 
            root.withdraw()
            info = tkinter.messagebox.showinfo("ATTENTION","This Regression model does not exist!!")
            root.destroy()
        #if both models exist the python evaluation is started
        else:
            fg = file_getter()  
            file_list = fg.get_automatic_list()
            
            remove_list = ['__header__', '__version__', '__globals__']
            drop_list = [
                            'GPS_time', 'q_T0', 'LaneChange_Approved',
                            'Left_Index', 'Right_Index', 'index',
                            'LatPos_abs', 'LongPos_abs', 'GPS_status',
                            'c02_left', 'c02_right', 'SteeringTorque',
                            'objectVelocityFrontRight', 'objectVelocityFrontLeft',
                            'objectAccelerationFrontRight', 'objectAccelerationFrontLeft',
                            'objectDistanceFrontRight', 'objectDistanceFrontLeft',
                        ]
            #the signals are imported according to the configuration from the signal_profiles json file
            fObj = open(py_eval_path+'/signal_profiles.json')
            signals = json.load(fObj)
            columns = signals[self.config+"_signals"]
            infer_file_list = []
            print("start conversion")
            for file in file_list:
                #the mat file is loaded from the _temp folder
                print(file)
                path = os.path.join(py_eval_path, "0_data_mat", file+".mat")
                df = pd.DataFrame()
                converter = conversion(path, df)
                mat_file = converter.read_mat()
                converter.save_converted_csv_to_temp(mat_file, file)
                #it is checked whether the signals from the signal profiles are contained in the mat file
                #skip file if it lacks signals 
                skip_file = converter.check_mat_signals(mat_file, columns)
                if skip_file:
                    print("INFO:",file, "will be skipped because the MAT file does not match the evaluation profile.")
                #convert to csv and save it if it has all the signals
                else:
                    #the button lane change marking is only done if the evaluation profile contains 'basic'  
                    if self.config == 'basic':
                        LC_raw = converter.laneChange_Cutter_raw()
                        df = converter.laneChange_Maker(LC_raw)
                    df = converter.mat_to_df(mat_file, columns, remove_list)
                    converter.save_to_csv(file)
                    infer_file_list.append(file)
            #infer all files that have sufficient signals
            print("start inference")
            for file in infer_file_list:
                print("evaluating", file)
                sample_rate = 10
                inference1 = inference(config=self.config)
                inference1.infer(class_model, reg_model, file, drop_list, sample_rate, signals=12)
                #start online evaluation comparing online predictions to kinetic labels and to regression prediction
                if self.config == 'online':
                    file_path = py_eval_path + '/1_full_sequence_dataset/' + file +'.csv'
                    df = pd.read_csv(file_path)
                    lc = df['laneChangeStateOnline']
                    directions = df['laneChangeDirectionOnline']
                    starts, ends = self.find_online_starts_ends(lc)
                    evaluator = automatic_evaluation(gentle_time_threshold=self.gentle_time_threshold, strict_time_threshold=self.gentle_time_threshold)
                    evaluator.evaluate_online(file, starts, ends, directions)
            #if the evaluation profile doesn't contain 'online' the regression prediction will be compared to the kinetic labels
            if not self.config == 'online':
                evaluator = automatic_evaluation(gentle_time_threshold=self.gentle_time_threshold,strict_time_threshold=self.strict_time_threshold)
                evaluator.evaluate(ground_truth='kinetic_labels', prediction='regression')
            #all the files generated during the python evaluation process will be deleted
            cleaner = cleanup()
            cleaner.post_inference_cleanup()
            print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')