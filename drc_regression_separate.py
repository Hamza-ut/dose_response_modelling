from drc.logistic import (
    fit_ll4,
    plot_covariance_heatmap,
    compute_parameter_errors,
    calculate_fit_metrics,
    calculate_confidence_intervals,
)
import pandas as pd

# DATA1
df = pd.read_csv("data/processed/ryegrass_metaremoved.csv")
x_data = df["conc"]
y_data = df["rootl"]

# # DATA2
# df = pd.read_csv("data/processed/lepidium_metaremoved.csv")
# x_data = df["conc"]
# y_data = df["weight"]
# four_parameter_logistic_regression(x_data, y_data)

# # DATA3
# df = pd.read_csv("data/processed/timepoint_vallo_9.5to10.5.csv")
# x_data = df["um"]
# y_data = df["od_normalized"]

# # DATA4
# df = pd.read_csv("data/processed/timepoint_sf_normalized.csv")
# x_data = df["xmic"]
# y_data = df["od_normalized"]


results = fit_ll4(x_data, y_data)
b_slope_fit, c_lower_fit, d_upper_fit, e_inflection_fit = results[0]
print(f"b(Slope)      = {b_slope_fit:.5f}")
print(f"c(Lower)      = {c_lower_fit:.5f}")
print(f"d(Upper)      = {d_upper_fit:.5f}")
print(f"e(Inflection) = {e_inflection_fit:.5f}")


std_errors, relative_std_errors = compute_parameter_errors(results[0], results[1])
print(
    f"Standard Errors(b,c,d,e): {std_errors[0]:.5f}, {std_errors[1]:.5f}, {std_errors[2]:.5f}, {std_errors[3]:.5f}"
)
print(
    f"Percentage/Relative Standard Errors(b,c,d,e): {relative_std_errors[0]:.2f}%, {relative_std_errors[1]:.2f}%, {relative_std_errors[2]:.2f}%, {relative_std_errors[3]:.2f}%"
)

plot_covariance_heatmap(results[1])

r_square, rmse = calculate_fit_metrics(y_data, x_data, results[0])
print(f"R²: {r_square:.5f}")
print(f"RMSE: {rmse:.5f}")


names = ["b(Slope)", "c(Lower)", "d(Upper)", "e(ED50/EC50)"]
confidence_percentage = 95
confidence_intervals = calculate_confidence_intervals(
    y_data, results[0], results[1], names, confidence_percentage
)
for interval in confidence_intervals:
    print(
        f"{interval['parameter']}: "
        f"{confidence_percentage}% Confidence Interval [{interval['lower_bound']:.5f} to {interval['upper_bound']:.5f}]"
    )
