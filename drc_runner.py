from drc.validation import (
    group_fields_validation,
    match_columns,
    validate_numeric_columns,
    validate_output_column,
)
from drc.io_utils import read_csv_file
from drc.correction import blank_correction, baseline_correction
from drc.normalization import normalization


def run_pipeline(
    file_path,
    group_fields,
    replicate_field,
    time_field,
    od_field,
    dose_field,
    well_field=None,  # None -> dataset has no well-identity column
    media_only_well_value=None,  # None -> no blank well -> use baseline_correction
    LOD=0.03,
    corrected_od_col_name="OD_corrected",
    normalized_od_col_name="OD_normalized",
):
    """Load, validate, correct (blank or baseline), and normalize one dataset.

    Which correction method runs is decided explicitly by whether
    `media_only_well_value` is given -- not auto-detected from the data.
    """
    # ---- 1. Validate + standardize the field-name arguments, once ----
    group_fields_validation(group_fields)

    group_fields = [f.strip().lower() for f in group_fields]
    replicate_field = replicate_field.strip().lower()
    time_field = time_field.strip().lower()
    od_field = od_field.strip().lower()
    dose_field = dose_field.strip().lower()
    well_field = well_field.strip().lower() if well_field is not None else None
    corrected_od_col_name = corrected_od_col_name.strip().lower()
    normalized_od_col_name = normalized_od_col_name.strip().lower()

    # ---- 2. Load (read_csv_file already lowercases the DataFrame's own columns) ----
    df = read_csv_file(file_path)

    # ---- 3. Validate the loaded data against those field names, once, upfront ----
    required_columns = group_fields + [replicate_field, dose_field, od_field, time_field]
    if well_field is not None:
        required_columns = required_columns + [well_field]

    match_columns(df, required_columns)
    validate_numeric_columns(df, [dose_field, od_field, time_field])
    validate_output_column(df, corrected_od_col_name)
    validate_output_column(df, normalized_od_col_name)

    # ---- 4. Pick the correction method explicitly, based on what was given ----
    if media_only_well_value is not None:
        corrected = blank_correction(
            df,
            group_fields=group_fields,
            replicate_field=replicate_field,
            time_field=time_field,
            od_field=od_field,
            dose_field=dose_field,
            media_only_well_value=media_only_well_value,
            LOD=LOD,
            corrected_od_col_name=corrected_od_col_name,
        )
    else:
        corrected = baseline_correction(
            df,
            group_fields=group_fields,
            replicate_field=replicate_field,
            time_field=time_field,
            od_field=od_field,
            dose_field=dose_field,
            well_field=well_field,
            LOD=LOD,
            corrected_od_col_name=corrected_od_col_name,
        )

    # ---- 5. Normalization is the same either way ----
    return normalization(
        corrected,
        group_fields=group_fields,
        replicate_field=replicate_field,
        time_field=time_field,
        od_field=corrected_od_col_name,  # normalize the CORRECTED column, not the raw one
        dose_field=dose_field,
        normalized_od_col_name=normalized_od_col_name,
    )


def main():
    # Dataset 1: has dedicated blank wells -> blank_correction
    out1 = run_pipeline(
        file_path="data/raw/timepoint_vallo.csv",
        group_fields=["Species"],
        replicate_field="Plt",
        time_field="Time_h",
        od_field="RawOD",
        dose_field="uM",
        media_only_well_value=-1,
        LOD=0.03,
    )
    out1.to_csv("data/processed/timepoint_vallo_normalized.csv", index=False)

    # Dataset 2: no blank wells, but has real well-level replication -> baseline_correction, well_field needed
    out2 = run_pipeline(
        file_path="data/raw/timepoint_sf.csv",
        group_fields=["Condition", "Ratio"],
        replicate_field="Plate",
        time_field="hour",
        od_field="Raw_od",
        dose_field="XMIC",
        well_field="Well",
        LOD=0.03,
    )
    out2.to_csv("data/processed/timepoint_sf_normalized.csv", index=False)

    # Dataset 3: no blank wells, and no well column at all -> baseline_correction, well_field omitted
    out3 = run_pipeline(
        file_path="data/raw/custom_growth_2-FMA_toxicity.csv",
        group_fields=["Species"],
        replicate_field="Replicate",
        time_field="Time",
        od_field="RawOD",
        dose_field="Dose",
        LOD=0.03,
    )
    out3.to_csv("data/processed/custom_growth_2-FMA_toxicity_normalized.csv", index=False)


if __name__ == "__main__":
    main()
