def drc_blank_correction(
    df,
    group_fields: list[str],
    time_field: str,
    od_field: str,
    dose_field: str,
    control_well_value: float | int = -1,
    LOD: float = 0.0,
    corrected_od_col_name: str = "OD_corrected",
):
    """Performs horizontal background subtraction using dedicated media blank wells.

    Parameters:
    -----------
    df : Input microplate assay dataset as pandas DataFrame.
    group_fields : list[str]
        Categorical grouping columns (e.g. ['Date', 'Plt']).
    time_field : str
        Column name representing time (e.g. 'Time_h').
    od_field : str
        Column name representing raw optical density readings.
    dose_field : str
        Column name representing compound concentration.
    control_well_value : float or int, default=-1
        Value in dose_field that designates sterile media blank wells.
    LOD : float, default=0.0
        Limit of Detection threshold floor.
    corrected_od_col_name : str, default='OD_corrected'
        Name of the newly created corrected column.

    Returns:
    --------
    pd.DataFrame: Copy of dataframe with background-subtracted optical density.
    """
    # 1. Type and emptiness checks
    if not isinstance(df, pd.DataFrame):
        raise TypeError(
            f"Expected a pandas DataFrame, but received {type(df).__name__}."
        )
    if df.empty:
        raise ValueError("Provided DataFrame is empty.")
    if not isinstance(group_fields, list):
        raise TypeError(
            f"`group_fields` must be a list of column names (e.g. ['Date']), "
            f"but received {type(group_fields).__name__}: {group_fields!r}."
        )
    if corrected_od_col_name in df.columns:
        raise ValueError(
            f"`corrected_od_col_name='{corrected_od_col_name}'` already exists as a "
            "column in the DataFrame. Choose a different name so it doesn't get "
            "silently overwritten."
        )

    # 2. Prevent mutating the caller's original DataFrame, and guarantee unique
    # row labels (e.g. pd.concat of multiple plates leaves duplicate labels,
    # which breaks the label-based assignment in the transform step below)
    df = df.copy()
    df = df.reset_index(drop=True)

    # 3. Construct grouping keys (deduplicating to prevent redundant keys)
    grouping_for_blanks = list(dict.fromkeys(group_fields + [time_field]))

    # 4. Check required columns (every column-name argument this function takes)
    required_cols = {*group_fields, time_field, od_field, dose_field}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise KeyError(
            f"Missing required column(s) in DataFrame: {sorted(missing_cols)}"
        )

    # 5. EDGE CASE GUARD: Check if control wells exist anywhere in dataset
    control_mask = df[dose_field] == control_well_value
    if not control_mask.any():
        raise ValueError(
            f"No blank wells found where `{dose_field} == {control_well_value}`. "
            "Check your dose_field or control_well_value parameter."
        )

    # 6. TRANSFORM LOGIC
    temp_col = "OD_minus1_of_each_timepoint"

    df[temp_col] = (
        df[df[dose_field] == control_well_value]
        .groupby(grouping_for_blanks)[od_field]
        .transform("mean")
    )

    df[temp_col] = df.groupby(grouping_for_blanks)[temp_col].transform("first")

    # 7. EDGE CASE GUARD: Check if any snapshot was missing a blank well
    if df[temp_col].isna().any():
        missing_count = df[temp_col].isna().sum()
        raise ValueError(
            f"{missing_count} rows belong to group snapshots that do not contain "
            f"a control well matching `{dose_field} == {control_well_value}`."
        )

    df[corrected_od_col_name] = df[od_field] - df[temp_col]
    df[corrected_od_col_name] = df[corrected_od_col_name].clip(lower=LOD)

    # 8. Clean up working column
    df.drop(columns=[temp_col], inplace=True)

    return df


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
    if not isinstance(group_fields, list):
        raise TypeError(
            f"`group_fields` must be a list of column names (e.g. ['Well']), "
            f"but received {type(group_fields).__name__}: {group_fields!r}."
        )
    if corrected_od_col_name in df.columns:
        raise ValueError(
            f"`corrected_od_col_name='{corrected_od_col_name}'` already exists as a "
            "column in the DataFrame. Choose a different name so it doesn't get "
            "silently overwritten."
        )

    # 2. Prevent mutating the caller's original DataFrame
    df = df.copy()

    # 3. Sort out grouping fields for baseline wells, ensuring uniqueness and preserving order
    grouping_for_baseline = list(dict.fromkeys(group_fields + [dose_field]))

    # 4. Check required columns
    required_cols = {*group_fields, time_field, od_field, dose_field}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise KeyError(
            f"Missing required column(s) in DataFrame: {sorted(missing_cols)}"
        )

    # 5. Chronological sort to ensure t=0 is the leading row for each well
    df = df.sort_values(by=grouping_for_baseline + [time_field])

    # 6. TRANSFORM LOGIC
    # NOTE: plain "first" silently skips NaN and would grab a later timepoint
    # as the baseline if t=0 is missing, so grab the literal first row instead.
    temp_col = "OD_zero_of_each_well"
    df[temp_col] = df.groupby(grouping_for_baseline)[od_field].transform(
        lambda well_od: well_od.iloc[0]
    )

    # 6b. EDGE CASE GUARD: Check if any well is missing a valid t=0 reading
    missing_mask = df[temp_col].isna()
    if missing_mask.any():
        missing_wells = (
            df.loc[missing_mask, grouping_for_baseline]
            .drop_duplicates()
            .to_dict("records")
        )
        raise ValueError(
            f"{len(missing_wells)} well(s) are missing a valid t=0 reading in "
            f"`{od_field}`, so no baseline can be computed for them: {missing_wells}"
        )

    # 7. Clip the corrected OD values to the Limit of Detection (LOD)
    df[corrected_od_col_name] = df[od_field] - df[temp_col]
    df[corrected_od_col_name] = df[corrected_od_col_name].clip(lower=LOD)

    # 8. Clean up working column
    df.drop(columns=[temp_col], inplace=True)

    return df


import pandas as pd

# Test for drc_baseline_correction function
df = pd.read_csv("data/raw/timepoint_sf.csv")
group_fields = ["Condition", "Ratio"]
dose_field = "XMIC"
od_field = "Raw_od"
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
new_data.to_csv("data/processed/timepoint_sf_baseline_corrected.csv", index=False)


# test for drc_blank_correction function
df = pd.read_csv("data/raw/timepoint_vallo.csv")
group_fields = ["Date", "Plt"]
dose_field = "uM"
od_field = "RawOD"
time_field = "Time_h"
# control well value and limit of detection of the instrument
control_well_value = -1
LOD = 0.03
new_data = drc_blank_correction(
    df,
    group_fields=group_fields,
    time_field=time_field,
    od_field=od_field,
    dose_field=dose_field,
    control_well_value=control_well_value,
    LOD=LOD,
    corrected_od_col_name="OD_corrected",
)
new_data.to_csv("data/processed/timepoint_vallo_blank_corrected.csv", index=False)
