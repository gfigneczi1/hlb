import scipy.io
import pandas as pd
import os
import glob
import pathlib
from matplotlib import pyplot as plt
import numpy as np

class conversion():

    def __init__(self, path, df = pd.DataFrame()):
        self.path = path
        self.df = df

    def read_mat(self):
        #path_to_mathfile = r'C:\Measurements\mat_files\for_eval\X590_20210407T080902Z_MP1307_DASy_0000.mat'
        mat_file = scipy.io.loadmat(self.path)
        return mat_file
    
    def check_mat_signals(self, mat_file, columns):
        skip_file = False
        keys = mat_file.keys()
        for name in columns:
            if name not in keys:
                skip_file = True
                print("WARNING:", name, "not in MAT file.")
                break
        return skip_file

    def mat_to_df(self, mat_file, columns, remove_list = []):
        keys = [key for key in columns]
        for key in keys:
            if key not in remove_list:
                self.df[key] = mat_file[key][0]
        self.df.dropna(axis='columns', inplace=True)
        self.df = self.df[columns]
        return self.df

    def save_converted_csv_to_temp(self, mat_file, file_name):
        df_temp = pd.DataFrame()
        keys = mat_file.keys()
        for key in keys:
            if key not in ["__header__","__version__","__globals__"]:
                df_temp[key] = mat_file[key][0]
        homedir = str(pathlib.Path(__file__).resolve().parent.parent.parent)
        file_path = os.path.join(homedir, "_temp", file_name+".csv")
        df_temp.to_csv(file_path)

    def data_plot(self, filename, save_path='plots'):
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        save_path = os.path.join(py_eval_path, save_path)
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        values = self.df.values
        # specify columns to plot
        groups = [i for i in range(self.df.shape[1])]
        #groups = [6, 17, 23, 9, 24]
        i = 1
        # plot each column
        f = plt.figure()
        f.set_figwidth(10)
        f.set_figheight(15)
        for group in groups:
            plt.subplot(len(groups), 1, i)
            plt.plot(values[:, group])
            plt.title(self.df.columns[group], y=0.5, loc='right')
            i += 1
        plt.savefig(os.path.join(save_path, filename+".pdf"))
        #plt.show()
        
    def save_to_csv(self, filename, save_path='1_full_sequence_dataset'):
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        save_path = os.path.join(py_eval_path, save_path)
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        self.df.to_csv(os.path.join(save_path, filename+".csv"), index=False)
        return
    
    # Cut the lane change section with respect to a certain time thresholds.
    def laneChange_Cutter_raw(self, threshold = 10000):
        temp = []
        flag_L = 0
        flag_R = 0
        start_index = 0
        for index in range(self.df.shape[0]):
            # if the blinker was set to on without any further instruction for a certain period of time, reset everything.
            if (flag_L == 1 or flag_R == 1) and index - start_index > threshold:
                flag_L = 0
                flag_R = 0
            
            # detect and notify the beginning of a lane change .   
            if self.df['Left_Index'][index] == 1 and flag_L == 0:
                if flag_R == 0:
                    start_index = index
                flag_L = 1
            elif self.df['Right_Index'][index] == 1 and flag_R == 0:
                if flag_L == 0:
                    start_index = index
                flag_R = 1
                
            # deals with interction where a right lane change is performed right after left lane change.
            if flag_L == 1 and flag_R == 1:
                
                
                # problem with this approach: continuous left-right turns will be mistakenly labelled.
                #end_index = index
                #period = self.df['GPS_time'][end_index] - self.df['GPS_time'][start_index]
                #if period < threshold:
                    #temp.append((start_index, end_index))
                
                flag_L = 0
                flag_R = 0
            
            # Mark the end of a lane change.
            if self.df['LaneChange_Approved'][index] == 1 and (flag_L == 1 or flag_R == 1):
                end_index = index
                period = self.df['GPS_time'][end_index] - self.df['GPS_time'][start_index]
                if period < threshold:
                    temp.append((start_index, end_index))
                flag_L = 0
                flag_R = 0
        return temp
    
    # Create a new column in data frame a labelling
    def laneChange_Maker(self, laneChange):
        temp = [0 for i in range(self.df.shape[0])]
        self.df['LaneChange'] = temp
        for start, end in laneChange:
            for i in range(start, end+1):
                self.df.at[i, 'LaneChange'] = 1
        return self.df
    
    # Resample the lane change section to make it more accurate using the yawRateESP.
    def yawRate_rasampling(df, LC, yaw_threshold, debounce):
        for start, end in LC:
            temp = 0
            count = 0
            for index in range(start, end+1):
                yawrate = df['yawRateESP'][index]
                if yawrate < yaw_threshold and count == 0:
                    temp = yawrate
                    count = 1
                    end = index
                    
                if count >= 1:
                    if abs(yawrate - temp) < 0.00001:
                        count+=1
                else:
                    count = 0
                
                if count*10 > debounce:
                    print("do something")
        return
    
    def merging_csv(self, files, filename):
        self.df = pd.concat(map(pd.read_csv, files), ignore_index=True)
        self.save_to_csv(filename, "dataset")