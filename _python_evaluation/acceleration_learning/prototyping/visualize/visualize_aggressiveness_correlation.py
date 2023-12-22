import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":
    curvatures_aggressiveness_df = pd.read_csv(r'C:\git\kdp_hlb_evalframework\_temp\analysis\curvature\scaled_curvature_summary.csv')
    curvatures_aggressiveness_df.sort_values(
        "driver_id",
        ascending=True
    )
    model_parameter_aggressiveness_df = pd.read_csv(r'C:\git\kdp_hlb_evalframework\_temp\analysis\model_parameters\scaled_parameters_summary.csv')
    model_parameter_aggressiveness_df.sort_values(
        "driver_id",
        ascending=True
    )
    plt.plot(curvatures_aggressiveness_df["curvature_mean"], model_parameter_aggressiveness_df["global_summary"], "bo")
    plt.title("Comparison of the curvature metric and the CA-PT2 model parameter profile per driver")
    plt.xlabel("Scaled aggressiveness from curvature metric")
    plt.ylabel("Scaled aggressiveness from model parameters profile")
    plt.show()
