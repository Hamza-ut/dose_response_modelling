def normalization(
    df,
    group_fields: list[str],
    replicate_field: str,
    time_field: str,
    od_field: str,
    dose_field: str,
    normalized_od_col_name: str = "OD_normalized",
):
    df = df.copy()

    # grouping keys for calculating mean baseline of untreated control wells (dose == 0)
    grouping_keys = group_fields + [replicate_field] + [time_field]

    # Mean baseline of untreated control wells (dose == 0)

    temp_col = "OD_zero_mean_of_each_group_replicate_timepoint"
    df[temp_col] = (
        df[df[dose_field] == 0].groupby(grouping_keys)[od_field].transform("mean")
    )

    # Broadcast baseline across each time snapshot
    df[temp_col] = df.groupby(grouping_keys)[temp_col].transform("first")

    # Divide and floor at 0.0
    df[normalized_od_col_name] = (df[od_field] / df[temp_col]).clip(lower=0.0)

    return df


# def main():
#     import pandas as pd

#     ## ---------------------------------------------------- ##
#     # DATA1
#     file_path = "data/processed/timepoint_vallo_blank_corrected.csv"
#     group_fields = ["Species"]
#     replicate_field = "Plt"
#     time_field = "Time_h"
#     od_field = "OD_corrected"
#     dose_field = "uM"
#     normalized_od_col_name = "OD_normalized"

#     df = pd.read_csv(file_path)

#     df = normalization(
#         df,
#         group_fields=group_fields,
#         replicate_field=replicate_field,
#         time_field=time_field,
#         od_field=od_field,
#         dose_field=dose_field,
#         normalized_od_col_name=normalized_od_col_name,
#     )
#     df.to_csv("data/processed/timepoint_vallo_normalized.csv", index=False)

#     ## ---------------------------------------------------- ##
#     # DATA2
#     file_path = "data/processed/timepoint_sf_baseline_corrected.csv"
#     group_fields = ["Condition", "Ratio"]
#     replicate_field = "Plate"
#     time_field = "hour"
#     od_field = "OD_corrected"
#     dose_field = "XMIC"
#     normalized_od_col_name = "OD_normalized"

#     df = pd.read_csv(file_path)

#     df = normalization(
#         df,
#         group_fields=group_fields,
#         replicate_field=replicate_field,
#         time_field=time_field,
#         od_field=od_field,
#         dose_field=dose_field,
#         normalized_od_col_name=normalized_od_col_name,
#     )
#     df.to_csv("data/processed/timepoint_sf_normalized.csv", index=False)

#     ## ---------------------------------------------------- ##
#     # DATA3
#     file_path = "data/processed/custom_growth_2-FMA_toxicity_baseline_corrected.csv"
#     group_fields = ["Species"]
#     replicate_field = "Replicate"
#     time_field = "Time"
#     od_field = "OD_corrected"
#     dose_field = "Dose"
#     normalized_od_col_name = "OD_normalized"

#     df = pd.read_csv(file_path)

#     df = normalization(
#         df,
#         group_fields=group_fields,
#         replicate_field=replicate_field,
#         time_field=time_field,
#         od_field=od_field,
#         dose_field=dose_field,
#         normalized_od_col_name=normalized_od_col_name,
#     )
#     df.to_csv("data/processed/custom_growth_2-FMA_toxicity_normalized.csv", index=False)


# if __name__ == "__main__":
#     main()
