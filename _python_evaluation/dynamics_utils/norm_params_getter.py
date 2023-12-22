import pandas as pd 
import numpy as np 
import os 
import pathlib

class norm_params_getter():
    def get_z_norm_params(self, metric):
        folderpath = str(pathlib.Path(__file__).resolve().parent)
        df_in = pd.read_csv(os.path.join(folderpath,"Norm_Parameters.csv"))
        for m, r, a, b, c, d in zip(df_in['metric'], df_in['Recording'], df_in['mu'], df_in['std'], df_in['min'], df_in['max']):
            if r == 'natural' and m == metric:
                mu = a
                if b == 0:
                    b = 0.0001
                std = b
                return mu, std
        return 0, 0
    def get_mm_norm_params(self, metric):
        folderpath = str(pathlib.Path(__file__).resolve().parent)
        df_in = pd.read_csv(os.path.join(folderpath,"Norm_Parameters.csv"))
        for m, r, a, b, c, d in zip(df_in['metric'], df_in['Recording'], df_in['mu'], df_in['std'], df_in['min'], df_in['max']):
            if r == 'natural' and m == metric: 
                if c == d:
                    d = c + 0.001
                n_min = c 
                n_max = d
                return n_min, n_max
        return 0, 0