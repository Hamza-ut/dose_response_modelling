def blank_correction(
    df,
    group_fields: list[str],
    replicate_field: str,
    time_field: str,
    od_field: str,
    dose_field: str,
    media_only_well_value: float | int = -1,
    LOD: float = 0.03,
    corrected_od_col_name: str = "OD_corrected",
):
    """Performs horizontal background subtraction using dedicated media blank wells.
    Parameters:
    df : Input microplate assay dataset as pandas DataFrame.
    group_fields : list[str]
        Categorical grouping columns (e.g. ['Date', 'Plt']).
    time_field : str
        Column name representing time (e.g. 'Time_h').
    od_field : str
        Column name representing raw optical density readings.
    dose_field : str
        Column name representing compound concentration.
    media_only_well_value : float or int, default=-1
        Value in dose_field that designates sterile media blank wells.
    LOD : float, default=0.03
        Limit of Detection threshold floor.
    corrected_od_col_name : str, default='OD_corrected'
        Name of the newly created corrected column.

    Returns:
    pd.DataFrame: Copy of dataframe with background-subtracted optical density.
    """
    # make a copy of the dataframe to avoid modifying the original
    df = df.copy()
    # Add replicate fields to the grouping
    grouping_for_blanks = group_fields + [replicate_field] + [time_field]

    # Check if control wells exist anywhere in dataset
    if not (df[dose_field] == media_only_well_value).any():
        raise ValueError(
            f"No control wells found in the dataset for `{dose_field} == {media_only_well_value}`."
        )

    # take mean of control well (-1) for each timepoint within each group and replicate
    temp_col = "OD_minus1_of_each_group_replicate_timepoint"

    df[temp_col] = (
        df[df[dose_field] == media_only_well_value]
        .groupby(grouping_for_blanks)[od_field]
        .transform("mean")
    )

    df[temp_col] = df.groupby(grouping_for_blanks)[temp_col].transform("first")

    # Make sure OD_minus1 shud have value for each row in the dataset
    if df[temp_col].isna().any():
        missing_count = df[temp_col].isna().sum()
        raise ValueError(
            f"{missing_count} rows belong to group snapshots that do not contain "
            f"a control well matching `{dose_field} == {media_only_well_value}`."
        )

    # subtract the OD of minus1 of each group(plate/time) to its corresponding OD value
    df[corrected_od_col_name] = df[od_field] - df[temp_col]

    # instrument limit of detection if any, otherwise default to 0.0
    df[corrected_od_col_name] = df[corrected_od_col_name].clip(lower=LOD)

    # return the corrected dataframe
    return df


def baseline_correction(
    df,
    group_fields: list[str],
    replicate_field: str,
    time_field: str,
    od_field: str,
    dose_field: str,
    well_field: str | None = None,
    LOD: float = 0.03,
    corrected_od_col_name: str = "OD_corrected",
):
    # make a copy of the dataframe to avoid modifying the original
    df = df.copy()

    # No well_field given -> the dataset has no well-identity column at all,
    # so group_fields + dose_field is the only grouping available (and is
    # assumed to already isolate one well). Only run the uniqueness check,
    # and only fall back to well-based grouping, when well_field is actually
    # given -- we can't check or group by a column that doesn't exist.
    if (
        well_field is not None
        and df.groupby(group_fields + [dose_field])[well_field].nunique().max() > 1
    ):
        grouping_for_baseline = group_fields + [replicate_field] + [well_field]
    else:
        grouping_for_baseline = group_fields + [dose_field]

    # take each well's own lowest OD reading across its whole time series
    temp_col = "OD_minimum_of_each_well"
    df[temp_col] = df.groupby(grouping_for_baseline)[od_field].transform("min")

    # EDGE CASE GUARD: a well with no valid OD readings at all has no minimum
    if df[temp_col].isna().any():
        missing_wells = (
            df.loc[df[temp_col].isna(), grouping_for_baseline]
            .drop_duplicates()
            .to_dict("records")
        )
        raise ValueError(
            f"{len(missing_wells)} well(s) have no valid `{od_field}` readings at all, "
            f"so no baseline can be computed for them: {missing_wells}"
        )

    # subtract each well's own minimum OD from all of its own readings
    df[corrected_od_col_name] = df[od_field] - df[temp_col]
    df[corrected_od_col_name] = df[corrected_od_col_name].clip(lower=LOD)

    return df


# def main():
#     import pandas as pd

#     ## ---------------------------------------------------- ##
#     # Blank correction
#     file_path = "data/raw/timepoint_vallo.csv"
#     group_fields = ["Species"]
#     replicate_field = "Plt"
#     well_field = "Well"
#     dose_field = "uM"
#     od_field = "RawOD"
#     time_field = "Time_h"
#     media_only_well_value = -1
#     LOD = 0.03  # because this insturment LOD was known
#     corrected_od_col_name = "OD_corrected"

#     df = pd.read_csv(file_path)

#     df = blank_correction(
#         df,
#         group_fields=group_fields,
#         replicate_field=replicate_field,
#         time_field=time_field,
#         od_field=od_field,
#         dose_field=dose_field,
#         media_only_well_value=media_only_well_value,
#         LOD=LOD,
#         corrected_od_col_name=corrected_od_col_name,
#     )
#     df.to_csv("data/processed/timepoint_vallo_blank_corrected.csv", index=False)

#     ## ---------------------------------------------------- ##
#     # Baseline correction on dataset with well information
#     file_path = "data/raw/timepoint_sf.csv"
#     group_fields = ["Condition", "Ratio"]
#     replicate_field = "Plate"
#     time_field = "hour"
#     od_field = "Raw_od"
#     dose_field = "XMIC"
#     well_field = "Well"
#     LOD = 0.03  # because this instrument LOD was not known so using default.
#     corrected_od_col_name = "OD_corrected"

#     df = pd.read_csv(file_path)

#     df = baseline_correction(
#         df,
#         group_fields=group_fields,
#         replicate_field=replicate_field,
#         well_field=well_field,
#         time_field=time_field,
#         od_field=od_field,
#         dose_field=dose_field,
#         LOD=LOD,
#         corrected_od_col_name=corrected_od_col_name,
#     )
#     df.to_csv("data/processed/timepoint_sf_baseline_corrected.csv", index=False)

#     ## ---------------------------------------------------- ##
#     # Baseline correction second dataset where there is no well information
#     file_path = "data/raw/custom_growth_2-FMA_toxicity.csv"
#     group_fields = ["Species"]
#     replicate_field = "Replicate"
#     time_field = "Time"
#     od_field = "RawOD"
#     dose_field = "Dose"
#     well_field = None  # no well information in this dataset
#     LOD = 0.03  # because this instrument LOD was not known so using default.
#     corrected_od_col_name = "OD_corrected"

#     df = pd.read_csv(file_path)

#     df = baseline_correction(
#         df,
#         group_fields=group_fields,
#         replicate_field=replicate_field,
#         well_field=well_field,
#         time_field=time_field,
#         od_field=od_field,
#         dose_field=dose_field,
#         LOD=LOD,
#         corrected_od_col_name=corrected_od_col_name,
#     )
#     df.to_csv(
#         "data/processed/custom_growth_2-FMA_toxicity_baseline_corrected.csv",
#         index=False,
#     )


# if __name__ == "__main__":
#     main()
