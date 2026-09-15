import numpy as np
import pandas as pd
from drc_baseline_correction import drc_baseline_correction
from drc_blank_correction import drc_blank_correction


from drc_four_pl_regression import four_parameter_logistic_regression


def fit_grouped_drc(df, group_cols, dose_col="XMIC", signal_col="OD_normalized"):
    """Fits 4PLL model independently for each unique combination of group columns."""
    results = []

    # Group by multiple columns (e.g., Condition and Ratio)
    for group_keys, group_df in df.groupby(group_cols):
        print(f"\n--- Fitting Curve for Group: {group_keys} ---")

        x_data = group_df[dose_col].values
        y_data = group_df[signal_col].values

        # Clean out any potential NaNs for the regression
        valid_idx = ~np.isnan(x_data) & ~np.isnan(y_data)
        x_data = x_data[valid_idx]
        y_data = y_data[valid_idx]

        if len(np.unique(x_data)) < 4:
            print(f"Skipping {group_keys}: Not enough unique dose points for 4PLL.")
            continue

        popt, pcov = four_parameter_logistic_regression(x_data, y_data)

        if popt is not None:
            b, c, d, e = popt
            # Handle group_keys format whether it's a tuple or single string
            if isinstance(group_cols, list) and len(group_cols) > 1:
                group_dict = dict(zip(group_cols, group_keys))
            else:
                group_dict = {
                    (
                        group_cols[0] if isinstance(group_cols, list) else group_cols
                    ): group_keys
                }

            results.append(
                {
                    **group_dict,
                    "b_slope": b,
                    "c_lower": c,
                    "d_upper": d,
                    "e_EC50": e,
                }
            )

    return pd.DataFrame(results)


# Dataset1
df = pd.read_csv("data/processed/timepoint_vallo_normalized.csv")

group_fields = ["Species"]
dose_field = "uM"
od_field = "OD_normalized"
time_field = "Time_h"

medium_only = -1
min_time = 9.5
max_time = 10.5


# # Dataset 2:
# df = pd.read_csv("data/processed/timepoint_sf_normalized.csv")

# group_fields = ["Condition", "Ratio"]
# dose_field = "XMIC"
# od_field = "OD_normalized"
# time_field = "hour"

# medium_only = -1
# min_time = 9.5
# max_time = 10.5


# Filter by your timepoint window and exclude medium only samples
filtered_df = df[
    (df[time_field] >= min_time)
    & (df[time_field] <= max_time)
    & (df[dose_field] != medium_only)
].copy()

summary_results = fit_grouped_drc(
    filtered_df, group_cols=group_fields, dose_col=dose_field, signal_col=od_field
)

print("\nFinal Summary Results:")
print(summary_results)
