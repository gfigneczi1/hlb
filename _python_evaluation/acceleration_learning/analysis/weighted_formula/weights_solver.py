import sys, os, json, pathlib
import logging
from collections import defaultdict
from typing import List
import matplotlib.pyplot as plt
import numpy as np
from acceleration_learning.preprocessing.situation_extraction.driving_situation_window import driving_situation_window
from acceleration_learning.preprocessing.situation_extraction.driving_event_enum import DrivingEventEnum
from acceleration_learning.exception.empty_window_exception import EmptyWindowsException
from acceleration_learning.analysis.weighted_formula.aggressiveness_formula import prepare_speedup_formula_parameters

class WeightsSolver:
    """Optimizer of the weights of the level of dynamic formula"""
    def __init__(self):
        pass

    def generate_features_dataset(
        self,
        speedup_windows: List[driving_situation_window],
        signal_profile: dict
    ) -> list:
        """Compute the features of speed up situations from the avt dataset and return them"""
        variables_per_speedup = []
        for speedup in speedup_windows:
            try:
                parameters = prepare_speedup_formula_parameters(speedup.data, signal_profile)
            except EmptyWindowsException:
                continue
            variables_per_speedup.append(parameters)
        return variables_per_speedup

    def __random_weights(self, inputs_dimension: int):
        return np.random.randn(inputs_dimension)

    def __compute_load(self, variables_list: list, weights: dict):
        loads = []
        for variables in variables_list:
            if len(variables) == len(weights):
                load = variables @ weights
                loads.append(load)
        return loads

    def __normalize_weights(self, weights: np.array):
        weights = weights - min(weights)
        return weights / sum(weights)

    def find_optimum(self, dataset_file: str):
        """Optimize the formula's weights to maximize the variance of the aggressiveness between 
        drivers"""
        variables_per_trip_dataset = dict()
        with open(dataset_file, "r", encoding='utf-8') as json_file:
            variables_per_trip_dataset = json.load(json_file)
        # variables_per_trip_dataset:
        # {
        #   trip1: [
        #       [j1a, p1a, ...],
        #       [j1b, p1b, ...],
        #       ...]
        # }
        n_iteration = 1000000
        mutation_rate = 0.0005
        standard_deviation_target = 1
        expected_accuracy = 0.01

        weights_len = len(next(iter(variables_per_trip_dataset.values()))[0])
        weights = []
        for _ in range(weights_len):
            weights.append(1/weights_len)
        weights = np.array(weights)
        logging.info("Initial weights: %s", weights)

        logs = {
            "error":[],
            "mean":[],
            "standard_deviation":[]
        }
        
        weights_log = defaultdict(list)

        variables_means_per_trip = []
        for trip_variables_list in variables_per_trip_dataset.values():
            variables_means = []
            variables_list = np.transpose(np.array(trip_variables_list))
            # variables_list: [
            #   [j1a, j1b, ...],
            #   [p1a, p1b, ...],
            #   ...]
            for variables in variables_list:
                variables_means.append(np.mean(variables))
            # variables_means: [
            #   j1_mean,
            #   p1_mean,
            #   ...]
            variables_means_per_trip.append(variables_means)

        for _ in range(n_iteration):
            # layer
            loads = self.__compute_load(variables_means_per_trip, weights)
            # error
            loads_standard_deviation = np.std(loads)
            loads_mean = np.mean(loads)

            error = standard_deviation_target - loads_standard_deviation  #+ loads_mean - mean_target
            error = np.sum(np.abs(error))
            if error < expected_accuracy:
                logging.info("Solution found")
                break
            mutated_weights = weights + (self.__random_weights(len(weights)) * mutation_rate)
            normalized_mutated_weights = self.__normalize_weights(mutated_weights)
            mutant_loads = self.__compute_load(variables_means_per_trip, normalized_mutated_weights)
            mutant_loads_standard_deviation = np.std(mutant_loads)
            #mutant_loads_mean = np.mean(mutant_loads)
            mutation_error = standard_deviation_target - mutant_loads_standard_deviation  #+ mutant_loads_mean - mean_target
            mutation_error = np.sum(np.abs(mutation_error))
            if mutation_error < error:
                weights = mutated_weights
            logs["error"].append(error)
            logs["mean"].append(loads_mean)
            logs["standard_deviation"].append(loads_standard_deviation)
            weights_log["Jerk"].append(weights[0])
            weights_log["pedal"].append(weights[1])
            weights_log["overshoot"].append(weights[2])
            weights_log["duration"].append(weights[3])
        
        for key, values in weights_log.items():
            plt.plot(values, label=key)
        plt.legend()
        plt.xlabel("Iteration")
        plt.ylabel(f"Weights")
        plt.title(f"Evolution of the weights during the training phase.")
        plt.savefig(r'C:\git\kdp_hlb_evalframework\_temp\analysis\weighted_formula\weights_evolution.png')
        plt.clf()   # Clear figure
        for key, values in logs.items():
            plt.plot(values, label=key)
        plt.legend()
        plt.xlabel("Iteration")
        plt.ylabel(f"Algorithm state")
        plt.title(f"Evolution of the tuning performances during the training phase. \nFinal error: {error:.2f}, final standard deviation: {loads_standard_deviation:.2f}")
        plt.savefig(r'C:\git\kdp_hlb_evalframework\_temp\analysis\weighted_formula\learning_rate.png')
        ##plt.show()
        
        logging.info(f"Error {error:.3f}, mean {loads_mean:.3f}, std {loads_standard_deviation:.3f}")
        weights = self.__normalize_weights(weights)
        logging.info("Optimized and normalized weights: %s", weights)
        return weights

    def save_weights(self, output_file_path: str, weights: list, situation: DrivingEventEnum):
        weights_labels = ["velocity_jerk_weight", "accelerator_pedal_position_weight", "velocity_overshoot_weight", 'duration_weight']
        speedup_dict = dict()
        for i, label in enumerate(weights_labels):
            speedup_dict[label] = weights[i]
        weights_dict = {situation.value: speedup_dict}
        with open(output_file_path, "w") as json_file:
            json.dump(weights_dict, json_file, indent=4)
