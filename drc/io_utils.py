import pandas as pd
from pathlib import Path


# CSV Reader
def read_csv_file(file_path: str | Path) -> pd.DataFrame:
    """
    Read a CSV file and create pandas df.

    :param file_path: Path to the CSV file
    :returns: Loaded DataFrame.
    """
    file_path = Path(file_path)
    # Check 1: file exists
    if not file_path.exists():
        raise FileNotFoundError(
            f"File '{file_path}' does not exist. Check your config file for the correct path."
        )
    # Check 2: only CSV allowed
    if file_path.suffix.lower() != ".csv":
        raise ValueError(
            f"Only .csv files are allowed: '{file_path.name}' not a CSV file."
        )
    # Try reading CSV with robust error handling
    try:
        df = pd.read_csv(file_path)
        df.columns = df.columns.str.strip().str.lower()
    except pd.errors.EmptyDataError:
        raise ValueError(f"CSV file '{file_path}' is empty.")
    except pd.errors.ParserError as e:
        raise ValueError(
            f"Failed to parse CSV file '{file_path}': {e}. "
            "Possible causes: bad delimiter, inconsistent columns, or malformed data."
        )
    except Exception as e:
        raise ValueError(f"CSV Reader Error for file '{file_path}': {e}")
    # Check 3: empty DataFrame after reading
    if df.empty:
        raise ValueError(f"CSV file '{file_path}' contains no data.")

    return df
