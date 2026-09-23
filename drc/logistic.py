from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.stats import t, linregress
import pandas as pd
import numpy as np


def ll4_formula(x, b, c, d, e):
    """4-Parameter Logistic Model formula.
    Parameters:
        x: Independent variable (dose/concentration)
        b: Slope (Hill slope)
        c: Lower limit (Bottom asymptote)
        d: Upper limit (Top asymptote)
        e: Inflection point(Effective dose / concentration)
    """
    return c + (d - c) / (1.0 + (x / e) ** b)


def fit_ll4(x_data, y_data, use_bounds=False):
    x_data = np.asarray(x_data, dtype=float)
    y_data = np.asarray(y_data, dtype=float)

    # read r_drc_vs_python_fitting.md to understand these initial guesses are chosen and why they are better then R_DRC in some cases
    # LOWER AND UPPER INITIAL
    c_lower_initial = np.percentile(y_data, 5)
    d_upper_initial = np.percentile(y_data, 95)

    # SLOPE INITIAL
    # get rid of doses of zero or medium wells otherwise log will throw error
    positive_mask = (
        (y_data > c_lower_initial) & (y_data < d_upper_initial) & (x_data > 0)
    )

    # EDGE CASE GUARD: too few points survive the trimming for a meaningful
    # regression (need at least 2 degrees of freedom: n_used - 2 parameters
    # being fit >= 2). Below that, linregress either returns a degenerate
    # zero-residual "perfect fit" (n_used == 2) or silently NaN (n_used < 2)
    # instead of raising -- catch it here with a clear message instead.
    n_used = positive_mask.sum()
    if n_used < 4:
        raise ValueError(
            f"Only {n_used} point(s) remain after excluding non-positive doses "
            f"and points outside the 5th-95th percentile range; need at least "
            f"4 for a meaningful slope regression. Dataset has {len(y_data)} "
            f"points total."
        )

    log_x = np.log(x_data[positive_mask])
    # subtract upper and lower limits from y_data and divide to turn them into ratio (compress the range), and finally take log of it to linearize it
    # this is basically inversing of 4PL formula
    transformed_y = np.log(
        (d_upper_initial - y_data[positive_mask])
        / (y_data[positive_mask] - c_lower_initial)
    )
    # take linear regression
    regression = linregress(log_x, transformed_y)
    b_slope_initial = regression.slope

    # INFLECTION INITIAL
    e_inflection_initial = np.exp(-regression.intercept / regression.slope)

    initial_guess = [
        b_slope_initial,
        c_lower_initial,
        d_upper_initial,
        e_inflection_initial,
    ]

    fit_kwargs = {"p0": initial_guess, "maxfev": 5000}
    if use_bounds:
        # setting use_bounds=True swaps scipy's algorithm from LM to TRF
        lower_bounds = [0, np.min(y_data) - 1, np.min(y_data), 0]
        upper_bounds = [100, np.max(y_data) + 1, np.max(y_data), np.max(x_data)]
        fit_kwargs["bounds"] = (lower_bounds, upper_bounds)

    popt, pcov = curve_fit(ll4_formula, x_data, y_data, **fit_kwargs)

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
