def drc_normalization(
    df,
    group_fields,
    time_field,
    od_field,
    dose_field,
    normalized_od_col_name: str = "OD_normalized",
):
    """Normalizes optical density by dividing treated well signal by the mean untreated control.

    Parameters:
    -----------
    df : pd.DataFrame
        Input dataset containing OD readings and dose information.
    group_fields : list[str]
        Categorical grouping columns defining unique plate snapshots (e.g., ['Date', 'Plt']).
    time_field : str
        Column name for the time axis.
    od_field : str
        Column name containing corrected OD readings to normalize.
    dose_field : str
        Column name containing drug concentration/dose.
    normalized_od_col_name : str, default='OD_normalized'
        Name of the resulting normalized column.

    Returns:
    --------
    pd.DataFrame: A copy of the DataFrame with the normalized column appended.
    """

    df = df.copy()

    # Deduplicate grouping keys
    grouping_keys = list(dict.fromkeys(group_fields + [time_field]))

    # 1. Mean baseline of untreated control wells (dose == 0)
    df["OD_zero_drug"] = (
        df[df[dose_field] == 0].groupby(grouping_keys)[od_field].transform("mean")
    )

    # 2. Broadcast baseline across each time snapshot
    df["OD_zero_drug"] = df.groupby(grouping_keys)["OD_zero_drug"].transform("first")

    # 3. Divide and floor at 0.0
    df[normalized_od_col_name] = (df[od_field] / df["OD_zero_drug"]).clip(lower=0.0)

    # Clean up working column
    df.drop(columns=["OD_zero_drug"], inplace=True)
    return df


import pandas as pd

df = pd.read_csv("timepoint_vallo_step111.csv")
group_fields = ["Date", "Plt"]
dose_field = "uM"
od_field = "OD_corrected"
time_field = "Time_h"


new_data = drc_normalization(
    df,
    group_fields=group_fields,
    time_field=time_field,
    od_field=od_field,
    dose_field=dose_field,
    normalized_od_col_name="OD_normalized",
)

new_data.to_csv("timepoint_vallo_normalized.csv", index=False)
