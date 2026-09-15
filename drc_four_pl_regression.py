from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.stats import t
import pandas as pd
import numpy as np


def ll4_formula(x, b, c, d, e):
    """4-Parameter Logistic Model formula.
    Parameters:
        x: Independent variable (dose/concentration)
        b: Slope (Hill slope)
        c: Lower limit (Bottom asymptote)
        d: Upper limit (Top asymptote)
        e: ED50 / EC50 (Effective dose / concentration)
    """
    return c + (d - c) / (1.0 + (x / e) ** b)


def fit_ll4(x_data, y_data):
    """Fits a 4-parameter logistic function and returns optimal parameters and covariance."""
    x_data = np.asarray(x_data, dtype=float)
    y_data = np.asarray(y_data, dtype=float)

    # Initial guesses mapped to (b, c, d, e)
    b_slope_initial = 1.0
    c_lower_initial = np.percentile(y_data, 5)
    d_upper_initial = np.percentile(y_data, 95)
    e_ed50_initial = np.median(x_data)
    initial_guess = [b_slope_initial, c_lower_initial, d_upper_initial, e_ed50_initial]

    # Bounds mapped to match [b, c, d, e] order
    # Lower bounds: [slope_min, lower_limit_min, upper_limit_min, ed50_min]
    # Upper bounds: [slope_max, lower_limit_max, upper_limit_max, ed50_max]
    lower_bounds = [0, np.min(y_data) - 1, np.min(y_data), 0]
    upper_bounds = [100, np.max(y_data) + 1, np.max(y_data), np.max(x_data)]
    bounds = (lower_bounds, upper_bounds)

    popt, pcov = curve_fit(
        ll4_formula,
        x_data,
        y_data,
        p0=initial_guess,
        bounds=bounds,
        maxfev=5000,
    )

    return popt, pcov


def plot_covariance_heatmap(cov_matrix):
    plt.figure(figsize=(7, 5))

    # 'coolwarm' maps negative numbers to blue, zero to white, and positive numbers to red
    plt.imshow(cov_matrix, cmap="coolwarm", interpolation="nearest")
    plt.colorbar(label="Raw Covariance")
    plt.xticks(
        range(4),
        ["b(Slope)", "c(Lower)", "d(Upper)", "e(ED50/EC50)"],
    )
    plt.yticks(
        range(4),
        ["b(Slope)", "c(Lower)", "d(Upper)", "e(ED50/EC50)"],
    )
    plt.title("Raw Covariance Matrix Heatmap")
    plt.show()


def compute_parameter_errors(popt, pcov):
    variance_values = np.diag(pcov)
    std_errors = np.sqrt(variance_values)
    relative_std_errors = std_errors / np.abs(popt) * 100

    return std_errors, relative_std_errors


def calculate_fit_metrics(y_true, x_data, popt):
    """Calculates R² and RMSE to evaluate overall 4PLL curve fit quality."""
    y_pred = ll4_formula(x_data, *popt)
    residuals = y_true - y_pred

    # R-squared calculation
    ss_res = np.sum(residuals**2)  # sum squared residuals
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r_square = 1.0 - (ss_res / ss_tot) if ss_tot != 0 else np.nan

    # RMSE calculation
    mse = np.mean(residuals**2)  # mean squared residuals
    rmse = np.sqrt(mse)

    return r_square, rmse


def calculate_confidence_intervals(
    y_data,
    popt,
    pcov,
    names=["b", "c", "d", "e"],
    confidence_percentage=95,
):
    """
    Calculates t-distribution confidence intervals.
    default parameter names are b,c,d,e and default confidence_percentage is 95
    """
    level = confidence_percentage / 100.0
    degrees_of_freedom = len(y_data) - len(popt)
    t_score = t.ppf((1 + level) / 2, df=degrees_of_freedom)

    standard_errors = np.sqrt(np.diag(pcov))

    intervals = []
    for name, est, error in zip(names, popt, standard_errors):
        intervals.append(
            {
                "parameter": name,
                "lower_bound": est - t_score * error,
                "upper_bound": est + t_score * error,
            }
        )

    return intervals


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
# x_data = df["uM"]
# y_data = df["OD_normalized"]


results = fit_ll4(x_data, y_data)
b_slope_fit, c_lower_fit, d_upper_fit, e_ed50_fit = results[0]
print(f"b(Slope)      = {b_slope_fit:.5f}")
print(f"c(Lower)      = {c_lower_fit:.5f}")
print(f"d(Upper)      = {d_upper_fit:.5f}")
print(f"e(ED50/EC50)  = {e_ed50_fit:.5f}")


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
