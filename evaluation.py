import matlab_evaluation
import os
from _python_evaluation.python_evaluation import python_evaluation
from _python_evaluation.train import Train

class Evaluate:

    def execute(config):
        if 'Python' in config['evaluation_profile']:
            if 'train' in config['evaluation_profile']:
                Train.execute(config['training_step'], config['classification_model'], config['regression_model'])
            elif 'online' in config['evaluation_profile']:
                py_eval = python_evaluation(config='online', strict_time_threshold=config['strict_threshold'], gentle_time_threshold=config['gentle_threshold'])
                py_eval.execute(config['classification_model'], config['regression_model'])
            elif 'light' in config['evaluation_profile']:
                py_eval = python_evaluation(config='light', strict_time_threshold=config['strict_threshold'], gentle_time_threshold=config['gentle_threshold'])
                py_eval.execute(config['classification_model'], config['regression_model'])
            else:
                py_eval = python_evaluation(config='basic', strict_time_threshold=config['strict_threshold'], gentle_time_threshold=config['gentle_threshold'])
                py_eval.execute(config['classification_model'], config['regression_model'])
            
        else:
            print('MATLAB evaluation started')
            # Working direction is set to path it over to matlab evaluations
            workingDirectory = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_temp")
            meas_eval = matlab_evaluation.initialize()
            meas_eval.matlab_evaluation_cover(workingDirectory, nargout=0)
            meas_eval.terminate()
