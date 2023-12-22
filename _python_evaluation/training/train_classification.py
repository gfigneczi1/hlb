import tensorflow as tf
from tensorflow.keras.optimizers import schedules, SGD, Adam
import keras
import datetime
import pandas as pd
import numpy as np
import os
import pathlib
from keras.preprocessing.image import ImageDataGenerator
import math
from scipy.signal import lfilter
from tensorflow.keras.utils import to_categorical
#from other directories
from _python_evaluation.models.classification_models import *
from _python_evaluation.utils.utilities import *
#global variables
initial_learning_rate = 0.001
RANDOM_SEED = 21
np.random.seed(RANDOM_SEED)
tf.random.set_seed(RANDOM_SEED)
'''
This class contains the classification
It is initialized with the number of epochs to train with, the batch size to train with and the input shape of the model
'''
class classification():
    
    def __init__(self, epochs, batch_size, input_shape, signals):
        self.epochs = epochs
        self.batch_size = batch_size
        self.input_shape = input_shape
        self.signals = signals
    '''
    This is a Learning rate schedular that decreases the learning rate every 30 epochs during training
    '''
    def scheduler(self, epoch, lr):
        drop_rate = 0.1
        epochs_drop = 30.0
        return initial_learning_rate * math.pow(drop_rate, math.floor(epoch/epochs_drop))
    '''
    This data loader returns the X data for the Classification training, validation and testing
    consisting of N files with 14 or less columns and 90 rows of data that is reshaped into shape (N, 36, 35, 1) 
    The y is simply a 1 if it is a lane change and a 0 if it is not
    '''
    def data_loading(self, data_path):
        X_data = []
        y_data = []
        util = utilities()
        for folder in os.listdir(data_path):
            sub_path = data_path + '/' + folder
            if folder != 'LeftRightTurn':
                for file in os.listdir(sub_path):
                    file_path = sub_path + '/' + file
                    df = pd.read_csv(file_path)
                    #df = util.data_filtering(df)
                    out = util.normalize(df)
                    if out.shape[0] < 90:
                        print(file_path)
                        print(df.shape)
                    out = out[:,0:self.signals]
                    X_data.append(out)
                    if folder == "LaneChange":
                        y_data.append(0)
                    elif folder == "Not":
                        y_data.append(1)
        data_X = np.array(X_data)
        data_y = np.array(y_data)
        data_X = data_X.reshape(data_X.shape[0], self.input_shape[0], self.input_shape[1], self.input_shape[2])
        print(data_X.shape)
        print(data_y.shape)
        data_y = to_categorical(data_y, num_classes=2)
        return data_X, data_y
    '''
    This method executes the training of the classification model. The model has to be input.
    '''
    def train(self, model, model_name):
        print("start classification training")
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        #the data loader is used to load the training, validation and testing data
        train_datagen = ImageDataGenerator(zoom_range=0.1, horizontal_flip=True)
        val_datagen = ImageDataGenerator()
        train_path=os.path.join(py_eval_path, "3_cut_dataset", "train")
        train_X, train_y = self.data_loading(train_path)
        train_it = train_datagen.flow(train_X, train_y, batch_size=self.batch_size)
        val_path=os.path.join(py_eval_path, "3_cut_dataset", "val")
        val_X, val_y = self.data_loading(val_path)
        val_it = val_datagen.flow(val_X, val_y, batch_size=self.batch_size)
        test_path=os.path.join(py_eval_path, "3_cut_dataset", "test")
        test_X, test_y = self.data_loading(test_path)
        #initialize optimizer
        sgd = SGD(learning_rate=initial_learning_rate, momentum=0.9)
        adam = Adam(learning_rate=initial_learning_rate)
        #compile model
        model.compile(optimizer=adam, loss='binary_crossentropy', metrics=[keras.metrics.Precision()])
        # 2 callbacks are used: the tensorboard callback for logging and the lr_callback to schedule the learning rate
        log_dir = py_eval_path+"/logs/fit_classification/" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir)
        early_callback = tf.keras.callbacks.EarlyStopping(
            monitor='val_accuracy', min_delta=0, patience=2, verbose=0,
            mode='auto', baseline=None, restore_best_weights=False
        )
        lr_callback = tf.keras.callbacks.LearningRateScheduler(self.scheduler, verbose=1)
        #fit model
        model.fit(train_it, steps_per_epoch=len(train_it),
            callbacks=[lr_callback, tensorboard_callback],
            validation_data=val_it, validation_steps=len(val_it), epochs=self.epochs, verbose=1)
        #make prediction yhat
        pre = []
        raw = []
        yhat = model.predict(test_X, verbose=0)
        #the predictions are two values between 0 and 1
        #if the first value is larger, the data is predicted to contain a lane change
        #otherwise the model is predicted to not contain a lane change
        for i in range(len(yhat)):
            if yhat[i][0] > yhat[i][1]:
                pre.append(1)
            else:
                pre.append(0)
            raw.append(yhat[i])
        #evaluation of classification prediction performance
        acc = 0
        for x, y, z in zip(pre, raw, test_y):
            print("Predicted result: {} with confident level {}. Ground truth: {}".format(x, y, z))
            if x == z[0]:
                acc+=1
        print('test accuracy: {:.2f}'.format(acc/len(test_y)))
        #the model is saved for later usage during inference
        model.save(os.path.join(py_eval_path, 'models', model_name+'.hdf5'))