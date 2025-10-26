# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
# wrapper_covert_csv_calc_lesion_percent.py

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Union, Dict, Any
import pandas as pd

# Import your two modules
from convert_csv_formats import convert_lesion_to_validation
from quantify_section_areas_and_volumes import (
    SeriesQuantConfig,
    quantify_section_areas_and_volumes,
    QuantifyOutput,
)

def run_convert_and_quantify(
    input_path: Union[str, Path],
    *,
    # Geometry/scaling knobs
    microns_per_pixel: float,                # REQUIRED (e.g., 568.0/500.0)
    section_thickness_um: float = 75.0,
    section_regex: str = r"sect(\d+)",

    # Converter knobs
    preserve_image_name: bool = True,        # keep "Image Name" so section numbers can be parsed
    pattern: str = "*_Lesion_Area.csv",
    suffix: str = "_converted_for_validation.csv",
    verbose_convert: bool = True,

    # Quantifier knobs
    plot_split_when_prompt: bool = True,     # show the striatum/lesion split plot before asking 1/2
    per_file_output_dir: Optional[Union[str, Path]] = None,  # if None, quantifier will prompt per file
    write_summary_json: bool = True,

    # Aggregation of tallies across files
    aggregate_output_dir: Optional[Union[str, Path]] = None,  # if provided and multiple files, write a summary CSV here
) -> Dict[str, Any]:
    """
    Pipeline:
      1) Convert input CSV(s) from Lesion_Area format to validation format (keeping 'Image Name' by default).
      2) Quantify areas/volumes for each converted CSV.

    Parameters
    ----------
    input_path : str | Path
        A single CSV file OR a directory containing matching CSVs (pattern controls matching).
    microns_per_pixel : float
        Spatial calibration in μm/px (e.g., if 500 px == 568 μm => 568.0/500.0).
    section_thickness_um : float
        Physical thickness per section (μm).
    section_regex : str
        Regex with a capturing group for section number from 'Image Name' (e.g., r"sect(\\d+)").
    preserve_image_name : bool
        Keep 'Image Name' during conversion (recommended so section numbers can be inferred).
    pattern : str
        When input_path is a directory, glob for files like '*_Lesion_Area.csv'.
    suffix : str
        Output suffix for converted files (directory mode).
    verbose_convert : bool
        Print per-file converter status.
    plot_split_when_prompt : bool
        Show split plot before prompting to choose inside/outside (1 or 2).
    per_file_output_dir : str | Path | None
        Where to write quantifier outputs. If None, the quantifier prompts per file (default behavior).
    write_summary_json : bool
        If True, quantifier writes a summary JSON per series.
    aggregate_output_dir : str | Path | None
        If provided and we process multiple files, write an aggregated tallies CSV here.

    Returns
    -------
    dict with keys:
        'converted_paths': List[Path]                   # converted CSVs
        'quantify_outputs': List[QuantifyOutput]        # outputs per converted CSV
        'aggregate_csv': Optional[Path]                 # path to aggregated tallies CSV (if written)
        'aggregate_df': Optional[pd.DataFrame]          # aggregated tallies dataframe (if created)
    """
    p = Path(input_path)

    # 1) Convert
    if p.is_dir():
        converted_paths = convert_lesion_to_validation(
            p,
            pattern=pattern,
            suffix=suffix,
            preserve_image_name=preserve_image_name,
            verbose=verbose_convert,
        )
    else:
        # For a single file, just let the converter decide the output name
        converted_paths = [
            convert_lesion_to_validation(
                p,
                output=None,  # None => "<stem>_converted_for_validation.csv"
                preserve_image_name=preserve_image_name,
            )
        ]

    # 2) Quantify each converted CSV
    quantify_outputs: List[QuantifyOutput] = []
    for conv_csv in converted_paths:
        cfg = SeriesQuantConfig(
            csv_path=conv_csv,
            microns_per_pixel=float(microns_per_pixel),
            section_thickness_um=float(section_thickness_um),
            section_regex=section_regex,
            plot_split_when_prompt=plot_split_when_prompt,
            output_dir=per_file_output_dir,   # if None, quantifier will prompt
            write_summary_json=write_summary_json,
        )
        q = quantify_section_areas_and_volumes(cfg)
        quantify_outputs.append(q)

    # 3) Aggregate tallies (optional, useful when processing many files)
    aggregate_df: Optional[pd.DataFrame] = None
    aggregate_csv: Optional[Path] = None

    if len(quantify_outputs) > 1:
        rows = []
        for q in quantify_outputs:
            row = dict(q.tallies)
            row["source_converted_csv"] = str(q.per_row_csv.name).replace("_areas_and_volumes.csv", "_converted_for_validation.csv")
            rows.append(row)
        aggregate_df = pd.DataFrame(rows)

        if aggregate_output_dir is not None:
            agg_dir = Path(aggregate_output_dir)
            agg_dir.mkdir(parents=True, exist_ok=True)
            aggregate_csv = agg_dir / "lesion_pipeline_summary.csv"
            aggregate_df.to_csv(aggregate_csv, index=False)
            print("Wrote aggregate summary:", aggregate_csv)

    return {
        "converted_paths": converted_paths,
        "quantify_outputs": quantify_outputs,
        "aggregate_csv": aggregate_csv,
        "aggregate_df": aggregate_df,
    }


# Optional: small CLI for quick terminal use
if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Convert Lesion_Area CSVs and quantify lesion volumes.")
    ap.add_argument("input", help="CSV file or directory containing '*_Lesion_Area.csv'.")
    ap.add_argument("--um-per-px", type=float, required=True, help="Microns per pixel (e.g., 568.0/500.0).")
    ap.add_argument("--thickness", type=float, default=75.0, help="Section thickness in μm (default 75).")
    ap.add_argument("--section-regex", default=r"sect(\d+)", help="Regex to extract section numbers from 'Image Name'.")
    ap.add_argument("--no-plot-split", action="store_true", help="Do NOT show the split plot before prompting.")
    ap.add_argument("--per-file-outdir", default=None, help="Output dir for quantifier. If omitted, will prompt per file.")
    ap.add_argument("--aggregate-outdir", default=None, help="If set and multiple files, write aggregate CSV here.")
    args = ap.parse_args()

    res = run_convert_and_quantify(
        args.input,
        microns_per_pixel=float(args.um_per_px),
        section_thickness_um=float(args.thickness),
        section_regex=args.section_regex,
        plot_split_when_prompt=(not args.no_plot_split),
        per_file_output_dir=args.per_file_outdir,
        aggregate_output_dir=args.aggregate_outdir,
    )
    print("Done. Converted:", len(res["converted_paths"]), "file(s).")

"""
from pathlib import Path
from wrapper_covert_csv_calc_lesion_percent import run_convert_and_quantify

# Calibration: if 500 px == 568 μm
um_per_px = 568.0 / 500.0

in_csv = Path("/Users/mirandahulsey-vincent/Desktop/compiled_lesion_ROIs/USA5510_upper_hemisphere_Lesion_Area.csv")

res = run_convert_and_quantify(
    in_csv,
    microns_per_pixel=um_per_px,
    section_thickness_um=50.0,
    section_regex=r"sect(\d+)",
    preserve_image_name=True,           # keep Image Name so section numbers can be parsed
    plot_split_when_prompt=True,        # show the split plot before the 1/2 prompt
    per_file_output_dir="/Users/mirandahulsey-vincent/Desktop/compiled_lesion_ROIs",
    aggregate_output_dir=None,          # not needed for single file
)

# Inspect the quantifier outputs
qout = res["quantify_outputs"][0]
print(qout.per_row_csv)
print(qout.tallies)


"""

