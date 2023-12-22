from _python_evaluation.training.train_preprocessing import Train_Preprocessing
from _python_evaluation.training.train_classification import classification
from _python_evaluation.training.train_regression import regression
from _python_evaluation.utils.cleanup import cleanup
from _python_evaluation.models.regression_models import simple_model
from _python_evaluation.models.classification_models import cifar_resnet20
from tkinter import * 
import tkinter.messagebox 
import pathlib
import os

class Train:
    def execute(step, classification_model_name, regression_model_name):
        print('Python training started')
        #The dataset needed to train the classification model and the Regression model will be generated here.
        if 'preprocessing' in step:
            print("=========PREPROCESSING===========")
            Train_Preprocessing.execute()                                
        #The classification model is trained in this step.
        elif 'classification' in step:
            #if the required dataset doesn't exist (preprocessing has not taken place).
            if os.path.exists(os.path.join(str(pathlib.Path(__file__).parent),'3_cut_dataset')):
                #check if a new model is trained or an old model has to be overwritten.
                model_path = os.path.join(str(pathlib.Path(__file__).parent), 'models', classification_model_name+'.hdf5')
                if os.path.isfile(model_path):
                    #send a warning message popup if an old model with the same name already exists.
                    root=Tk() 
                    root.withdraw()
                    result=tkinter.messagebox.askquestion('ATTENTION!','Do you want to overwrite %s?' % classification_model_name)
                    #if yes is selected in the popup, a new model will be trained and overwrite the old model
                    if result=='yes':
                        print("====CLASSIFICATION TRAINING======")
                        input_shape = (36, 30, 1) 
                        initial_learning_rate = 0.001
                        model=cifar_resnet20()
                        #the specifications of how the model is trained can be changed here
                        classification_trainer = classification(epochs=100, batch_size=16, input_shape=input_shape, signals=12)
                        classification_trainer.train(model, classification_model_name)
                    else:
                        info = tkinter.messagebox.showinfo("INFO","The Classification model training was canceled")
                    root.destroy()
                #no old model exists so a new model will be trained and saved
                else:
                    print("====CLASSIFICATION TRAINING======")
                    input_shape = (36, 30, 1) 
                    initial_learning_rate = 0.001
                    model=cifar_resnet20()
                    #the specifications of how the model is trained can be changed here
                    classification_trainer = classification(epochs=100, batch_size=16, input_shape=input_shape, signals=12)
                    classification_trainer.train(model, classification_model_name)
            else:
                root=Tk() 
                root.withdraw()
                info = tkinter.messagebox.showinfo("ATTENTION","The dataset needed to train the Classification model does not exist!! Generate it by executing the preprocessing step.")
                root.destroy()
        #The classification model is trained in this step.
        elif 'regression' in step:
            #if the required dataset doesn't exist (preprocessing has not taken place).
            if os.path.exists(os.path.join(str(pathlib.Path(__file__).parent),'4_relabeled_dataset')):
                #check if a new model is trained or an old model has to be overwritten.
                model_path = os.path.join(str(pathlib.Path(__file__).parent), 'models', regression_model_name+'.hdf5')
                if os.path.isfile(model_path):
                    #send a warning message popup if an old model with the same name already exists.
                    root=Tk() 
                    root.withdraw()
                    result=tkinter.messagebox.askquestion('ATTENTION!','Do you want to overwrite %s?' % regression_model_name)
                    #if yes is selected in the popup, a new model will be trained and overwrite the old model
                    if result=='yes':
                        print("=====REGRESSION TRAINING========")
                        input_shape = (90, 1, 2)
                        model=simple_model(input_shape)
                        #the specifications of how the model is trained can be changed here
                        regressor = regression(epochs=80, batch_size=16, input_shape=input_shape)
                        regressor.train(model, regression_model_name)
                    else:
                        info = tkinter.messagebox.showinfo("INFO","The Regression model training was canceled")
                    root.destroy()
                else:
                    print("=====REGRESSION TRAINING========")
                    input_shape = (90, 1, 2)
                    model=simple_model(input_shape)
                    #the specifications of how the model is trained can be changed here
                    regressor = regression(epochs=80, batch_size=16, input_shape=input_shape)
                    regressor.train(model, regression_model_name)
            #no old model exists so a new model will be trained and saved
            else:
                root=Tk() 
                root.withdraw()
                info = tkinter.messagebox.showinfo("ATTENTION","The dataset needed to train the Regression model does not exist!! Generate it by executing the preprocessing step.")
                root.destroy()
        #in this step the datasets generated for the training of the models (in the _python_evaluation folder) are all deleted.
        #this will fail if something is manually deleted
        #the logs from the training process will also be deleted
        elif 'clean' in step:
            #warning message popup where one has to confirm if one really wants to delete the datasets
            root=Tk() 
            root.withdraw()
            result=tkinter.messagebox.askquestion('ATTENTION!','Do you want to delete the training datasets?')
            if result=='yes':
                cleaner = cleanup()
                cleaner.post_training_cleanup()
                print("cleanup succeeded")
            root.destroy()
