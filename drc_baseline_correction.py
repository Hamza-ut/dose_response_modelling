def drc_baseline_correction(
    df,
    group_fields: list[str],
    time_field: str,
    od_field: str,
    dose_field: str,
    LOD: float = 0.0,
    corrected_od_col_name: str = "OD_corrected",
):
    """Performs vertical baseline subtraction down each well's time series using t=0 readings.

    Parameters:
    -----------
    df : pd.DataFrame
        Input microplate assay dataset.
    group_fields : list[str]
        Categorical grouping columns defining unique wells (e.g. ['Species', 'Replicate']).
    time_field : str
        Column name representing time (e.g. 'Time').
    od_field : str
        Column name representing raw optical density readings.
    dose_field : str
        Column name representing compound concentration.
    LOD : float, default=0.0
        Limit of Detection threshold floor.
    corrected_od_col_name : str, default='OD_corrected'
        Name of the newly created corrected column.

    Returns:
    --------
    pd.DataFrame: Copy of dataframe with baseline-subtracted optical density.
    """
    # 1. Type and emptiness checks
    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            f"Expected a pandas DataFrame, but received {type(df).__name__}."
        )
    if df.empty:
        raise ValueError("Provided DataFrame is empty.")

    # 2. Prevent mutating the caller's original DataFrame
    df = df.copy()

    # 3. Sort out grouping fields for baseline wells, ensuring uniqueness and preserving order
    grouping_for_baseline = list(dict.fromkeys(group_fields + [dose_field]))

    # 4. Check required columns
    required_cols = set(grouping_for_baseline + [od_field, dose_field])
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise KeyError(
            f"Missing required column(s) in DataFrame: {sorted(missing_cols)}"
        )

    # 5. Chronological sort to ensure t=0 is the leading row for each well
    df = df.sort_values(by=grouping_for_baseline + [time_field])

    # 6. TRANSFORM LOGIC
    temp_col = "OD_zero_of_each_well"

    df[temp_col] = df.groupby(grouping_for_baseline)[od_field].transform("first")

    df[corrected_od_col_name] = df[od_field] - df[temp_col]
    df[corrected_od_col_name] = df[corrected_od_col_name].clip(lower=LOD)

    # 8. Clean up working column
    df.drop(columns=[temp_col], inplace=True)

    return df


import pandas as pd

df = pd.read_csv("custom_growth_2-FMA_toxicity.csv")
group_fields = ["Species", "Replicate"]
dose_field = "Dose"
od_field = "RawOD"
time_field = "Time"

# limit of detection of the instrument
LOD = 0.0

new_data = drc_baseline_correction(
    df,
    group_fields=group_fields,
    time_field=time_field,
    od_field=od_field,
    dose_field=dose_field,
    LOD=LOD,
    corrected_od_col_name="OD_corrected",
)

new_data.to_csv("custom_baseline_corrected.csv", index=False)
