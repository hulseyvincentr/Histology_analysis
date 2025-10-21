# -*- coding: utf-8 -*-
# convert_csv_formats.py
from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional, Union, Tuple, overload
import re
import pandas as pd

__all__ = [
    "TARGET_COLS",
    "convert_dataframe",
    "convert_file",
    "convert_dir",
    "convert_lesion_to_validation",
]

"""
Convert CSVs from the 'USA5510_Lower_Hemisphere_Lesion_Area' format
into the 'sample_shapes_for_validation' format—usable as a library function.

- Drops the 'Image Name' column by default (you can keep it via preserve_image_name=True).
- Ensures exact target column order.
- Parses coordinate strings like "1818, 1818, ..., 1830," (trailing commas/spaces OK).
- If a region is marked as 'none present' (or equivalent), writes 0.0 to match your sample.

Target column order:
    Striatum_x_coords, Striatum_y_coords,
    Lesion_x_coords,   Lesion_y_coords,
    Area_X_x_coords,   Area_X_y_coords,
    LMAN_x_coords,     LMAN_y_coords
"""

# Target column order to match the example "sample_shapes_for_validation.csv"
TARGET_COLS: List[str] = [
    "Striatum_x_coords", "Striatum_y_coords",
    "Lesion_x_coords",   "Lesion_y_coords",
    "Area_X_x_coords",   "Area_X_y_coords",
    "LMAN_x_coords",     "LMAN_y_coords",
]

# Tokens interpreted as "no shape present" (compare after lowercasing)
_NONE_TOKENS = {"none", "none present", "nan", "", "none_present", "no shape", "not present", "n/a", "na"}


def _parse_coords_cell(cell: Union[str, float, int]) -> Optional[List[float]]:
    """
    Accepts a cell that is either:
      - a string with comma-separated numbers (possibly with spaces / trailing comma / parentheses), or
      - one of the "none" tokens above.

    Returns a list of floats, or None if no shape is present.
    """
    if cell is None:
        return None

    s = str(cell).strip()
    if s.lower() in _NONE_TOKENS:
        return None

    # Remove surrounding parentheses/brackets to normalize "(1, 2, 3)"-style strings.
    s = s.strip("()[]")

    # Split on commas; tolerate trailing commas and extra spaces.
    parts = [p.strip() for p in s.split(",") if p.strip() != ""]
    nums: List[float] = []

    for p in parts:
        # Extract the first number-like token in each piece (handles "123", "123.45", "-1.2e3")
        m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", p)
        if m:
            try:
                nums.append(float(m.group(0)))
            except Exception:
                # Ignore pieces that don't parse cleanly
                pass

    return nums if nums else None


def _format_as_tuple_string(vals: Optional[List[float]]) -> Union[str, float]:
    """
    For present shapes: return a tuple-like string "(v1, v2, v3, ...)".
    For missing shapes: return 0.0 (to mirror the sample file behavior).
    """
    if not vals:
        return 0.0

    clean = [v for v in vals if v is not None]

    def _fmt(v: float) -> str:
        as_int = int(v)
        return str(as_int) if abs(v - as_int) < 1e-9 else str(v)

    return "(" + ", ".join(_fmt(v) for v in clean) + ")"


def convert_dataframe(df: pd.DataFrame, *, preserve_image_name: bool = False) -> pd.DataFrame:
    """
    Convert a DataFrame in the USA5510 'Lesion_Area' wide format to the
    'sample_shapes_for_validation' style.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with columns like:
          ['Image Name', 'Striatum_x_coords', 'Striatum_y_coords', 'Lesion_x_coords',
           'Lesion_y_coords', 'Area_X_x_coords', 'Area_X_y_coords', 'LMAN_x_coords', 'LMAN_y_coords']
    preserve_image_name : bool
        If True and 'Image Name' is present, keep that column as the first column in the output.

    Returns
    -------
    pd.DataFrame
        Converted DataFrame with columns in TARGET_COLS order (and optionally 'Image Name' first).
    """
    # Normalize exact column names (strip whitespace)
    col_map = {col: col.strip() for col in df.columns}
    df = df.rename(columns=col_map)

    # Validate required columns
    missing = [c for c in TARGET_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Input is missing expected columns: {missing}")

    out_rows: List[dict] = []
    for _, row in df.iterrows():
        out_row = {}

        # Optionally preserve 'Image Name' as a leading column
        if preserve_image_name and "Image Name" in df.columns:
            out_row["Image Name"] = row["Image Name"]

        # Convert target columns in exact order
        for c in TARGET_COLS:
            parsed = _parse_coords_cell(row[c])
            out_row[c] = _format_as_tuple_string(parsed)

        out_rows.append(out_row)

    # Determine final column order
    if preserve_image_name and "Image Name" in df.columns:
        cols = ["Image Name"] + TARGET_COLS
    else:
        cols = TARGET_COLS

    return pd.DataFrame(out_rows, columns=cols)


def convert_file(in_path: Union[str, Path],
                 out_path: Optional[Union[str, Path]] = None,
                 *,
                 preserve_image_name: bool = False) -> Path:
    """
    Convert a single CSV file and write it out.

    Returns the output path.
    """
    in_path = Path(in_path)
    if out_path is None:
        out_path = in_path.with_name(in_path.stem + "_converted_for_validation.csv")
    out_path = Path(out_path)

    df = pd.read_csv(in_path)
    out_df = convert_dataframe(df, preserve_image_name=preserve_image_name)
    out_df.to_csv(out_path, index=False)
    return out_path


def convert_dir(in_dir: Union[str, Path],
                pattern: str = "*_Lesion_Area.csv",
                out_suffix: str = "_converted_for_validation.csv",
                *,
                preserve_image_name: bool = False,
                verbose: bool = True) -> List[Path]:
    """
    Batch-convert all matching CSV files in a directory.

    Returns a list of output paths that were successfully written.
    """
    in_dir = Path(in_dir)
    outputs: List[Path] = []
    for p in sorted(in_dir.glob(pattern)):
        out_p = p.with_name(p.stem + out_suffix)
        try:
            convert_file(p, out_p, preserve_image_name=preserve_image_name)
            outputs.append(out_p)
            if verbose:
                print(f"[OK] {p.name} -> {out_p.name}")
        except Exception as e:
            if verbose:
                print(f"[WARN] Skipping {p.name}: {e}")
    return outputs


# ──────────────────────────────────────────────────────────────────────────────
# Unified function-first API
# ──────────────────────────────────────────────────────────────────────────────

@overload
def convert_lesion_to_validation(input_path: Union[str, Path],
                                 *,
                                 output: Optional[Union[str, Path]] = ...,
                                 preserve_image_name: bool = ...,
                                 pattern: str = ...,
                                 suffix: str = ...,
                                 verbose: bool = ...
                                 ) -> Path: ...
@overload
def convert_lesion_to_validation(input_path: Union[str, Path],
                                 *,
                                 output: None = ...,
                                 preserve_image_name: bool = ...,
                                 pattern: str = ...,
                                 suffix: str = ...,
                                 verbose: bool = ...
                                 ) -> Union[Path, List[Path]]: ...

def convert_lesion_to_validation(input_path: Union[str, Path],
                                 *,
                                 output: Optional[Union[str, Path]] = None,
                                 preserve_image_name: bool = False,
                                 pattern: str = "*_Lesion_Area.csv",
                                 suffix: str = "_converted_for_validation.csv",
                                 verbose: bool = True
                                 ) -> Union[Path, List[Path]]:
    """
    Unified, callable entry point.

    If `input_path` is a FILE:
        - Writes to `output` if provided, otherwise uses "<stem>_converted_for_validation.csv".
        - Returns Path to the output CSV.

    If `input_path` is a DIRECTORY:
        - Globs by `pattern`, writes each with `<stem><suffix>.csv`.
        - Returns List[Path] of outputs.

    Parameters
    ----------
    input_path : str | Path
        A single CSV file OR a directory to batch-convert.
    output : str | Path | None
        Only used when `input_path` is a single file.
    preserve_image_name : bool
        Keep 'Image Name' column first if present.
    pattern : str
        Glob pattern used when `input_path` is a directory.
    suffix : str
        Output suffix appended to each stem when `input_path` is a directory.
    verbose : bool
        Print per-file status in directory mode.

    Returns
    -------
    Path | List[Path]
        Output path(s) written.
    """
    p = Path(input_path)
    if p.is_dir():
        return convert_dir(p, pattern=pattern, out_suffix=suffix,
                           preserve_image_name=preserve_image_name, verbose=verbose)
    else:
        return convert_file(p, output, preserve_image_name=preserve_image_name)
 

# ──────────────────────────────────────────────────────────────────────────────
# Optional: tiny CLI wrapper stays for convenience, but the function above is primary.
# ──────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Convert Lesion_Area CSVs to sample_shapes_for_validation format.")
    ap.add_argument("input", help="Input CSV file OR a directory to batch-convert.")
    ap.add_argument("-o", "--output", default=None,
                    help="Output CSV path (only used when input is a single file).")
    ap.add_argument("--pattern", default="*_Lesion_Area.csv",
                    help="Glob pattern when input is a directory.")
    ap.add_argument("--suffix", default="_converted_for_validation.csv",
                    help="Output suffix when input is a directory.")
    ap.add_argument("--preserve-image-name", action="store_true",
                    help="Keep 'Image Name' column as the first column if present.")
    ap.add_argument("--quiet", action="store_true", help="Suppress per-file messages in directory mode.")
    args = ap.parse_args()

    outs = convert_lesion_to_validation(
        args.input,
        output=args.output,
        preserve_image_name=args.preserve_image_name,
        pattern=args.pattern,
        suffix=args.suffix,
        verbose=not args.quiet,
    )
    if isinstance(outs, list):
        print(f"Converted {len(outs)} file(s).")
        for op in outs: print(" -", op)
    else:
        print("Wrote:", outs)


"""
# Preferred: import the function directly (works if the folder is on sys.path)
from convert_csv_formats import convert_lesion_to_validation
from pathlib import Path
import pandas as pd

# --- Single file ---
in_csv  = Path("/Volumes/my_own_ssd/microscope_images/USA5510/lower_hemisphere_ROIs/USA5510_Lower_Hemisphere_Lesion_Area.csv")
out_csv = in_csv.with_name(in_csv.stem + "_converted_for_validation.csv")

res = convert_lesion_to_validation(
    in_csv,
    output=out_csv,            # optional; if omitted, suffix is auto-added
    preserve_image_name=False, # set True to keep "Image Name" as first column
)
print("Wrote:", res)
display(pd.read_csv(res).head())

# --- Batch (directory) ---
in_dir = Path("/Volumes/my_own_ssd/microscope_images/USA5510/lower_hemisphere_ROIs")
outs = convert_lesion_to_validation(
    in_dir,
    pattern="*_Lesion_Area.csv",
    suffix="_converted_for_validation.csv",
    preserve_image_name=False,
    verbose=True,
)
print("Converted:", len(outs), "files")
for p in outs: print(" -", p)



"""