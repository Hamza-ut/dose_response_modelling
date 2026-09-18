import pandas as pd


# Fields validator
def group_fields_validation(group_fields: list[str]):
    if not isinstance(group_fields, list):
        raise TypeError(
            f"`group_fields` must be a list of column names, even if its a single column (e.g. ['Species']),"
            f"but received {type(group_fields).__name__}: {group_fields!r}."
        )


# Columns matcher
def match_columns(df: pd.DataFrame, required_columns: list[str]) -> None:
    """
    Ensure required columns exist in the DataFrame.

    :param df: Input DataFrame
    :param required_columns: Columns that must exist
    """
    # standardize column names for case-insensitive matching (without mutating the caller's DataFrame)
    df_columns = [col.strip().lower() for col in df.columns]
    required_columns = [col.strip().lower() for col in required_columns]

    missing = [col for col in required_columns if col not in df_columns]
    if missing:
        raise ValueError(
            "CSV file is missing required column(s) as specified in the config file. \n"
            f"Missing required column(s) from the CSV: {missing}. \n"
        )


# Columns validator
def validate_numeric_columns(df: pd.DataFrame, numeric_columns: list[str]) -> None:
    """
    Validate numeric columns for non-numeric or missing data.

    :param df: Input DataFrame
    :param numeric_columns: Columns expected to contain numeric data
    """
    # standardize column names for case-insensitive matching (without mutating the caller's DataFrame)
    column_lookup = {col.strip().lower(): col for col in df.columns}
    numeric_columns = [col.strip().lower() for col in numeric_columns]

    for col in numeric_columns:
        # Check 1: Existence
        if col not in column_lookup:
            raise KeyError(f"Column '{col}' not found in DataFrame.")

        series = df[column_lookup[col]]

        # Check 2: Numeric Type
        if not pd.api.types.is_numeric_dtype(series):
            raise TypeError(f"Column '{col}' must be numeric.")

        # Check 3: Missing Values
        if series.isna().any():
            raise ValueError(f"Column '{col}' contains NaN values. Clean data first.")

        # Check 4: Variance (The 'Identical' check)
        if series.nunique() <= 1:
            raise ValueError(f"Column '{col}' has identical values; cannot process.")


# Output column validator
def validate_output_column(df: pd.DataFrame, output_col_name: str) -> None:
    """
    Ensure the output column name doesn't already exist in the DataFrame

    :param df: Input DataFrame
    :param output_col_name: Column name about to be written into
    """
    if output_col_name in df.columns:
        raise ValueError(
            f"`{output_col_name}` already exists as a column in the DataFrame. "
            "Choose a different name so it doesn't get silently overwritten."
        )
