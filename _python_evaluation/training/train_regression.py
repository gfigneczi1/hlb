import tensorflow as tf
import numpy as np
import pandas as pd
import math
import os
import pathlib
import datetime
from _python_evaluation.models.regression_models import *
from _python_evaluation.utils.utilities import *
#global variables
initial_learning_rate = 0.001
RANDOM_SEED = 21
np.random.seed(RANDOM_SEED)
tf.random.set_seed(RANDOM_SEED)
'''
This class contains the regression
It is initialized with the number of epochs to train with, the batch size to train with and the input shape of the model
'''
class regression():

    def __init__(self, epochs, batch_size, input_shape=(90, 1, 2)):
        self.epochs = epochs
        self.batch_size = batch_size
        self.input_shape = input_shape
    '''
    This is a Learning rate schedular that decreases the learning rate by *0.1 every 20 epochs during training
    '''
    def scheduler(self, epoch, lr):
        drop_rate = 0.1
        epochs_drop = 20.0
        return initial_learning_rate * math.pow(drop_rate, math.floor(epoch/epochs_drop))
    '''
    This data loader returns the X-data for the Regression training, validation and testing
    consisting of N files with 2 columns and 90 rows of data that is reshaped into shape (N, 90, 1, 2) 
    The X-data consists of the relative lateral acceleration and the average c01
    The ground truth is the final kinetic labelling generated during the relabelling step
    '''
    def data_loading(self, data_path):
        X_data = []
        y_data = []
        util = utilities()
        for file in os.listdir(data_path):
            file_path = data_path + '/' + file
            df = pd.read_csv(file_path)

            accY = util.get_accY(df)
            c01 = util.get_c01_average(df)
            out = np.column_stack((accY, c01))
            X_data.append(out)

            lc3 = df['kinetic_LaneChange']
            start, end = util.find_start_end(lc3)
            y_data.append((start, end))
        data_X = np.array(X_data)
        data_y = np.array(y_data)
        #reshape: (N, 90, 1, 2) with N samples, 90 in the time dimension, 1 data column per channel with 2 channels
        data_X = data_X.reshape(data_X.shape[0], self.input_shape[0], self.input_shape[1], self.input_shape[2])
        print(data_X.shape)
        print(data_y.shape)
        return data_X, data_y
    '''
    This method performs the Regression training. The model has to be input.
    '''
    def train(self, model, model_name):
        print("start regression training")
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        #data for training, validation and testing is loaded
        train_path = os.path.join(py_eval_path,"4_relabeled_dataset","train","LaneChange")
        val_path = os.path.join(py_eval_path,"4_relabeled_dataset","val","LaneChange")
        test_path = os.path.join(py_eval_path,"4_relabeled_dataset","test","LaneChange")

        train_X, train_y = self.data_loading(train_path)
        val_X, val_y = self.data_loading(val_path)
        test_X, test_y = self.data_loading(test_path)
        #compile the model
        model.compile(optimizer='adam', loss='mse')
        # 2 callbacks are used: the tensorboard callback for logging and the lr_callback to schedule the learning rate
        lr_callback = tf.keras.callbacks.LearningRateScheduler(self.scheduler, verbose=1)
        log_dir = py_eval_path+"/logs/fit_regression/" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir)
        #fit model
        model.fit(train_X, train_y, epochs=self.epochs, batch_size=self.batch_size, validation_data=(val_X, val_y),
                        shuffle=True, callbacks=[lr_callback, tensorboard_callback])
        #prediction
        yhat = model.predict(test_X)
            
        total_start = 0
        total_end = 0
        total_mse_start = 0
        total_mse_end = 0
        big_start_errors = 0
        big_end_errors = 0
        ground_truth_negative = 0
        for i in range(yhat.shape[0]):
            print('----------------------------')
            print('prediction ', yhat[i,:], yhat[i,1]-yhat[i,0])
            print('ground_truth ', test_y[i,:], test_y[i,1]-test_y[i,0])
            
            length_difference = (yhat[i,1]-yhat[i,0])-(test_y[i,1]-test_y[i,0])
            print('differences ', yhat[i,:] - test_y[i,:], length_difference)
            if test_y[i, 0] == 0 and test_y[i,1] == 90:
                ground_truth_negative += 1
            else:
                if abs(yhat[i,0] - test_y[i,0]) >= 3.0:
                    big_start_errors += 1
                if abs(yhat[i,1] - test_y[i,1]) >= 3.0:
                    big_end_errors += 1 
            total_start += 10*abs(yhat[i,0] - test_y[i,0])
            total_end += 10*abs(yhat[i,1] - test_y[i,1])
            total_mse_start += (test_y[i,0] - yhat[i,0])*(test_y[i,0] - yhat[i,0])
            total_mse_end += (test_y[i,1] - yhat[i,1])*(test_y[i,1] - yhat[i,1])
        print('----------------------------')
        print('test mean absolute error: ', total_start/yhat.shape[0], total_end/yhat.shape[0])
        print('test mean squared error: ', total_mse_start/yhat.shape[0], total_mse_end/yhat.shape[0])
        print('number of files:', yhat.shape[0])
        print('ground truth negatives: ', ground_truth_negative)
        print('start off by more than 300ms:', big_start_errors, 'times')
        print('end off by more than 300ms:', big_end_errors, 'times')
        positives = yhat.shape[0] - ground_truth_negative
        print('accuracy of start:', (positives-big_start_errors)/positives)
        print('accuracy of end:', (positives-big_end_errors)/positives)
        model.save(os.path.join(py_eval_path,'models', model_name+'.hdf5'))