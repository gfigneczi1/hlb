"""Python module to get the overshoot of a speedup"""
import pandas as pds

def measure_overshoot(speed_up_df: pd.DataFrame) -> float:
    """Return the overshoot of the velocity"""
    initial_velocity = speed_up_df["LongitudinalVelocity_0x344_20"].iloc[0]
    final_velocity = speed_up_df["LongitudinalVelocity_0x344_20"].iloc[-1]
    velocity_change = final_velocity - initial_velocity
    maximal_velocity = max(speed_up_df["LongitudinalVelocity_0x344_20"])
    velocity_overshoot = maximal_velocity - final_velocity
    velocity_overshoot_ratio = velocity_overshoot / velocity_change
    return velocity_overshoot_ratio
