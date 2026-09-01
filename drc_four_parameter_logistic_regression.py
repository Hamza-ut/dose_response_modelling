import pandas as pd
import numpy as np
from scipy.optimize import curve_fit


def four_pl_formula(x, A, B, C, D):
    return A + (B - A) / (1.0 + (x / C) ** D)


# manual 4 Paramter Log-Logistic Function for DRC Modeling
def four_paramter_logistic_regression(x_data, y_data):
    """Fit a 4-parameter logistic function to the data.
    Input:
        x_data: array-like, independent variable (e.g., concentration)
        y_data: array-like, dependent variable (e.g., response)
    """
    # Convert inputs to clean numpy arrays
    x = np.asarray(x_data, dtype=float)
    y = np.asarray(y_data, dtype=float)

    # 1. Filter out NaN values and non-positive doses (x <= 0 causes zero-division in SciPy)
    mask = (x > 0) & (~np.isnan(x)) & (~np.isnan(y))
    x_fit = x[mask]
    y_fit = y[mask]

    if len(x_fit) < 4:
        raise ValueError(
            "At least 4 positive dose points are required to fit a 4PL model."
        )

    # 2. Initial guesses matching your percentile strategy
    A_guess = np.percentile(y_fit, 5)  # Bottom asymptote
    B_guess = np.percentile(y_fit, 95)  # Top asymptote
    C_guess = np.median(x_fit)  # Midpoint / EC50 guess
    D_guess = 1.0  # Hill slope

    initial_guess = [A_guess, B_guess, C_guess, D_guess]

    # 3. # Define bounds for parameters
    bounds = (
        [np.min(y_data) - 1, np.min(y_data), 0, 0],
        [np.max(y_data) + 1, np.max(y_data), np.max(x_data), 100],
    )
    # # bounds = ([0, 0, 0, 0], [np.inf, np.inf, np.inf, 10])

    try:
        popt, pcov = curve_fit(
            four_pl_formula,
            x_fit,
            y_fit,
            p0=initial_guess,
            bounds=bounds,
            maxfev=5000,
        )
        A_fit, B_fit, C_fit, D_fit = popt

        # R-squared quality calculation
        residuals = y_fit - four_pl_formula(x_fit, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else np.nan

        print(
            f"Fitted Parameters:\n  Bottom (A) = {A_fit:.6f}\n  Top (B)    = {B_fit:.6f}\n  EC50 (C)   = {C_fit:.6f}\n  Slope (D)  = {D_fit:.6f}\n  R²         = {r_squared:.6f}"
        )
        return popt, pcov

    except Exception as e:
        print(f"Curve fitting failed: {e}")
        return None, None


# df = pd.read_csv("timepoint_vallo_normalized.csv")

# # filtering for Time_h between 9.5 and 10.5 h gives a desired result.
# filtered_df = df[(df["Time_h"] > 9.5) & (df["Time_h"] < 10.5)]
# filtered_df = filtered_df[filtered_df["uM"] != -1]
# filtered_df = filtered_df[["uM", "Fit", "OD_normalized"]]


# x_data = filtered_df["uM"]
# y_data = filtered_df["OD_normalized"]

df = pd.read_csv("data/processed/ryegrass_metaremoved.csv")

x_data = df["rootl"]
y_data = df["conc"]

four_paramter_logistic_regression(x_data, y_data)
