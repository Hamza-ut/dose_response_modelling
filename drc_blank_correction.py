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

    # 2. Prevent mutating the caller's original DataFrame
    df = df.copy()

    # 3. Construct grouping keys (deduplicating to prevent redundant keys)
    grouping_for_blanks = list(dict.fromkeys(group_fields + [time_field]))

    # 4. Check required columns
    required_cols = set(grouping_for_blanks + [od_field, dose_field])
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


import pandas as pd

df = pd.read_csv("timepoint_vallo.csv")
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

new_data.to_csv("timepoint_vallo_step111.csv", index=False)
