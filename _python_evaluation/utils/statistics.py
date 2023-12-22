import matplotlib.pyplot as plt
import pandas as pd

def plot_hist_df(df: pd.DataFrame, column: str):
    """Display the data frame of the given column as a distribution"""
    plt.hist(df[column], bins="auto")
    plt.ylabel("Number of samples")
    plt.xlabel(column)
    plt.title("Distribution of samples of " + column)
    plt.show()

def plot_hist(data: list, name):
    """Display the data list's distribution"""
    plt.hist(data, bins="auto")
    plt.ylabel(f"Number of {name}")
    plt.xlabel(name)
    plt.title(f"Distribution of samples of {name}")
    plt.show()

def filter_outliers(df: pd.DataFrame, column: str, removal_percentage: float):
    """Remove the given percentage of the extremes values from the mean of the given column"""
    quantile_low = df[column].quantile(removal_percentage)
    quantile_high  = df[column].quantile(1-removal_percentage)
    df = df[(df[column] < quantile_high) & (df[column] > quantile_low)]  
    return df

def filter_outliers_all_columns(df: pd.DataFrame, removal_percentage: float):
    """Remove the given percentage of the extremes values from the mean of all columns"""
    for column in df:
        filter_outliers(df, column, removal_percentage)
    return df
