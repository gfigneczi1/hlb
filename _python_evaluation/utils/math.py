import pandas as pd

def value_is_in_tolerance(value, tolerance, offset=0) -> bool:
    """Check if the given value is (bounded) in the tolerance  around the offset"""
    return abs(value) <= offset + tolerance

def min_max_normalize(value: float, min_value: float, max_value: float) -> float:
    """Return the value normalized with the min max normalization
    (ZeroDivision exception is raised if the max value is zero)"""
    return (value - min_value) / max_value

def min_max_normalize_data_frame(data_frame: pd.DataFrame, column_to_normalize: str) -> pd.DataFrame:
    """Return the min max normalized data frame
    (ZeroDivision exception is raised if the max value is zero)"""
    scaled_data_frame = data_frame.copy()
    for column_name in data_frame.columns:
        if column_name == column_to_normalize:
            min_value = min(data_frame[column_to_normalize])
            max_value = max(data_frame[column_to_normalize])
            scaled_data_frame[column_to_normalize] = (data_frame[column_to_normalize] - min_value) / (max_value - min_value)
    return scaled_data_frame
