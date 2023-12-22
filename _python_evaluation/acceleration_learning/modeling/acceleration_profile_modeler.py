import os
import json
import logging
import numpy as np
from scipy.optimize import fmin
from scipy import signal
import pandas as pd
from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.modeling.model import Model
from acceleration_learning.modeling.model_parameter import CAPT2ModelParameter

class acceleration_profile_modeler():
    def trip_parameter_estimation(self, input_windows_directory_path, signal_profile, model="CA_PT2"):
        '''returns the parameters averaged over all speedups in a trip and the parameters for the individual speedups in the trip
        
        keyword arguments:
        input_windows_directory_path -- path to the folder containing the extracted speedup data for the trip
        '''
        all_speedup_parameters = []
        speedup_index = 0
        for file in os.listdir(input_windows_directory_path):
            if "speed" in file:
                normalized_longitudinal_velocity = self.get_velocity_data(os.path.join(input_windows_directory_path, "speed_up{index}.csv".format(index = speedup_index)), signal_profile)
                parameters_i, include_in_average = self.speedup_parameter_estimation(normalized_longitudinal_velocity, model)
                if include_in_average:
                    all_speedup_parameters.append(parameters_i)
                speedup_index += 1
        all_speedup_parameters = np.array(all_speedup_parameters)
        if len(all_speedup_parameters) > 0:
            trip_parameters = np.mean(all_speedup_parameters, axis=0)
        else:
            trip_parameters = np.array([np.nan])
        return trip_parameters, all_speedup_parameters

    def get_velocity_data(self, speedup_data_path, signal_profile):
        '''returns acceleration signal from one speedup normalized between 0 and 1 and sampled down to 100 samples
        
        keyword arguments:
        speedup_data_path -- path to the speedup that needs to be analyzed
        '''
        AVT_data_table = pd.read_csv(speedup_data_path)
        longitudinal_velocity = list(AVT_data_table[signal_profile["Longitudinal_Velocity"]].values)
        first_second = min(30, len(longitudinal_velocity))
        velocity_change_difference = max(longitudinal_velocity) - min(longitudinal_velocity[1:first_second])
        normalized_longitudinal_velocity = (longitudinal_velocity - min(longitudinal_velocity[1:first_second]))/ velocity_change_difference
        normalized_longitudinal_velocity = signal.resample(normalized_longitudinal_velocity, 100)
        return normalized_longitudinal_velocity

    def G_PT1(self, x):
        '''returns the transfer function of the P2U model as a scipy.signal.lti
        the transfer function is: G(s) = 1 / (Tp1*s + 1)

        keyword argument:
        x -- the parameters x0 and D formatted like this: x = [Tp1]
        '''
        num = [1]
        den = [x[0], 1]
        return signal.lti(num, den)

    def G_PT2(self, x):
        '''returns the transfer function of the P2U model as a scipy.signal.lti
        the transfer function is: G(s) = 1 / ((1/w0**2)*s**2 + (2*D/w0)*s + 1)

        keyword argument:
        x -- the parameters x0 and D formatted like this: x = [w0, D]
        '''
        num = [1]
        den = [1/x[0]**2, 2*x[1]/x[0], 1]
        return signal.lti(num, den)

    def get_ca_pt2_element(self, x, U):
        '''returns the ca_pt2 estimated y_CA_PT2(t=U)

        keyword argument:
        x -- the parameters x0 and D formatted like this: x = [D, w0, tc, alfa, t_delay]
        U -- time point for which the CA_PT2 estimated function is taken
        '''
        # returns output value based on time = U
        D = x[0]
        T = 1/x[1]
        tc = x[2]
        alfa = x[3]
        tdelay = x[4]
        v0 = alfa*tc
        G = 1

        if U < tdelay:
            Y = 0
        elif U < (tc+tdelay):
            Y = alfa*(U-tdelay)
        else:
            t = U-tc-tdelay
            if D > 1:
                l1 = (-D+np.sqrt(D**2-1))/T
                l2 = (-D-np.sqrt(D**2-1))/T
                c1 = (v0-G)*(l2-alfa)/(l2-l1)
                c2 = (alfa-l1*(v0-G))/(l2-l1)
                Y = c1*np.exp(l1*t) + c2*np.exp(l2*t) + G
            elif D == 1:
                l = -1/T
                c1 = v0-G
                c2 = alfa - l*(v0-G)
                Y = c1*np.exp(l*t) + c2*t*np.exp(l*t) + G
            else:
                A = -D/T
                B = np.sqrt(1-D**2)/T
                c1 = v0-G
                c2 = (alfa-A*(v0-G))/B
                Y = np.exp(A*t)*(c1*np.cos(B*t) + c2*np.sin(B*t)) + G

        return Y

    def speedup_parameter_estimation(self, normalized_longitudinal_velocity, model="CA_PT2"):
        '''
        returns - the estimated parameters for one individual speedup
                    - PT1 model: [Tp1]
                    - PT2 model: [w0, D]
                    - CA_PT2 model: [D, w0, tc, alfa, t_delay]
                - boolean include_in_average that is False for exception cases where the optimization did not converge.
        Here the speedup is modeled as the step response of the P2U model G(x). The optimal parameters of w0 and D are found using fmin s.t. the mean squared error
        between the estimated signal and the normalized longitudinal velocity of the speedup is minimized.
        If the optimization method does not converge, the NaN will be returned

        keyword arguments: 
        normalized_longitudinal_velocity -- velocity signal from one speedup normalized between 0 and 1 and sampled down to 100 samples
        model -- the name of the model to estimate with
        '''
        y = np.array(normalized_longitudinal_velocity)
        t = np.arange(100)
        def PT1_fun(x):
            '''returns the mean squared error between the measurement velocity signal and the signal estimated with the PT1 model

            keyword argument:
            x -- the parameters x0 and D formatted like this: x = [Tp1]
            '''
            _, yx = signal.step(self.G_PT1(x), T=t)
            return ((yx - y)**2).mean()

        def PT2_fun(x):
            '''returns the mean squared error between the measurement velocity signal and the signal estimated with the PT2 model

            keyword argument:
            x -- the parameters x0 and D formatted like this: x = [w0, D]
            '''
            _, yx = signal.step(self.G_PT2(x), T=t)
            return ((yx - y)**2).mean()

        def CA_PT2_fun(x):
            '''returns the mean squared error between the measurement velocity signal and the signal estimated with the CA_PT2 model

            keyword argument:
            x -- the parameters x0 and D formatted like this: x = [D, w0, tc, alfa, t_delay]
            '''
            yx = []
            for time_point in t:
                Y = self.get_ca_pt2_element(x, time_point)
                yx.append(Y)
            return ((yx - y)**2).mean()

        try:
            if model == "PT1":
                # Set initial guess for the parameters [Tp1]
                x0 = [3.5]
                # Use fmin to find the minimum
                res = fmin(PT1_fun, x0)
            elif model == "PT2":
                # Set initial guess for the parameters [w0, D]
                x0 = [0.5,0.7]
                # Use fmin to find the minimum
                res = fmin(PT2_fun, x0)
            else:
                # Set initial guess for the parameters [D, w0, tc, alfa, t_delay]
                x0 = [0.5, 1, 3, 0.125, 0.5]
                # Use fmin to find the minimum
                res = fmin(CA_PT2_fun, x0)
                logging.info(f"The estimated parameters [D, w0, tc, alfa, t_delay] are: {res}")
            return res, True
        except:
            #return initial guess if it does not converge
            return np.nan, False

    def save_trip_parameter_output(
        self,
        output_parameter_file_path: str,
        trip_parameters, model=Model.CONSTANT_ACCELERATION_AND_SECOND_ORDER_STEP_RESPONSE
    ):
        '''saves the average parameters of a trip to a json file

        keyword arguments:
        output_parameter_file_path -- the name of a trip
        trip_parameters -- the averages of all the parameters in a trip
        model -- the name of the model used for estimation
        '''
        if model == Model.FIRST_ORDER_STEP_RESPONSE:
            parameter_dictionary = {"Tp1": trip_parameters[0]}
        elif model == Model.SECOND_ORDER_STEP_RESPONSE:
            parameter_dictionary = {"w0": trip_parameters[0], "D": trip_parameters[1]}
        else:
            parameter_dictionary = {
                CAPT2ModelParameter.D.value: trip_parameters[0],
                CAPT2ModelParameter.W.value: trip_parameters[1],
                CAPT2ModelParameter.TC.value: trip_parameters[2],
                CAPT2ModelParameter.A.value: trip_parameters[3],
                CAPT2ModelParameter.TD.value: trip_parameters[4]
            }
        with open(output_parameter_file_path, 'w', encoding="utf-8") as new_file:
            json.dump(parameter_dictionary, new_file, indent=4)
