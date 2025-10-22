# -*- coding: utf-8 -*-
# quantify_section_areas_and_volumes.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union

import os
import re
import json
import pandas as pd
import matplotlib.pyplot as plt

# Requires shapely: pip install shapely
try:
    from shapely.geometry import Polygon, LineString
    from shapely.ops import split as shapely_split
except Exception as e:
    raise ImportError("This module requires 'shapely'. Install with: pip install shapely") from e


# ──────────────────────────────────────────────────────────────────────────────
# Helpers: parsing & geometry
# ──────────────────────────────────────────────────────────────────────────────
_NONE_TOKENS = {
    "none", "none present", "none_present", "no shape", "not present",
    "na", "n/a", "nan", "", "0", "0.0"
}

def _parse_coords_cell(cell: Union[str, float, int, None]) -> Optional[List[float]]:
    """Parse '(1,2,3)', '1, 2, 3,', or numbers; return list[float] or None for 'none' tokens."""
    if cell is None:
        return None
    if isinstance(cell, (int, float)):
        if float(cell) == 0.0:
            return None
        return [float(cell)]
    s = str(cell).strip()
    if s.lower() in _NONE_TOKENS:
        return None
    s = s.strip("()[]")
    parts = [p.strip() for p in s.split(",") if p.strip()]
    out: List[float] = []
    for p in parts:
        m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", p)
        if m:
            try:
                out.append(float(m.group(0)))
            except Exception:
                pass
    return out or None


def _coords_to_polygon(x: Optional[List[float]], y: Optional[List[float]]) -> Optional[Polygon]:
    if not x or not y:
        return None
    n = min(len(x), len(y))
    if n < 3:
        return None
    return Polygon(zip(x[:n], y[:n]))


def _coords_to_line(x: Optional[List[float]], y: Optional[List[float]]) -> Optional[LineString]:
    if not x or not y:
        return None
    n = min(len(x), len(y))
    if n < 2:
        return None
    return LineString(zip(x[:n], y[:n]))


def _get_section_number(image_name: str, pattern: str) -> Optional[int]:
    m = re.search(pattern, str(image_name))
    return int(m.group(1)) if m else None


# ──────────────────────────────────────────────────────────────────────────────
# Public config, outputs, and API
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class SeriesQuantConfig:
    csv_path: Union[str, Path]
    microns_per_pixel: float            # e.g., if 500 px == 568 μm, set 568/500
    section_thickness_um: float = 75.0  # thickness per section in μm
    section_regex: str = r"sect(\d+)"   # how to extract the section number
    plot_split_when_prompt: bool = True # show plot before asking 1/2
    output_dir: Optional[Union[str, Path]] = None  # ask user if None
    write_summary_json: bool = True     # also write totals JSON


@dataclass
class QuantifyOutput:
    dataframe: pd.DataFrame
    per_row_csv: Path
    tallies: Dict[str, float]
    summary_json: Optional[Path]


def _choose_output_dir(cfg: SeriesQuantConfig) -> Path:
    """Resolve the output directory by config or user prompt, then create it."""
    if cfg.output_dir is None:
        default_dir = Path(cfg.csv_path).parent
        user = input(f"Enter output directory (press Enter for default: {default_dir}): ").strip()
        out_dir = Path(os.path.expanduser(user)) if user else default_dir
    else:
        out_dir = Path(os.path.expanduser(str(cfg.output_dir)))
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir.resolve()


def _tally_totals(df: pd.DataFrame, results_csv_path: Path) -> Dict[str, float]:
    """Compute final totals/percentages, matching your keys."""
    def colsum(col: str) -> float:
        return float(df[col].sum()) if col in df.columns else 0.0

    total_area_x_volume = colsum("Area X Volume (um^3)")
    total_overlap_volume_x = colsum("Overlap Volume with Area X (um^3)")
    total_lman_volume = colsum("LMAN Volume (um^3)")
    total_overlap_volume_lman = colsum("Overlap Volume with LMAN (um^3)")
    total_lesion_outside_striatum = colsum("Lesion Volume outside Striatum (um^3)")

    percent_area_x_lesioned = (100.0 * total_overlap_volume_x / total_area_x_volume) if total_area_x_volume else 0.0
    percent_lman_lesioned = (100.0 * total_overlap_volume_lman / total_lman_volume) if total_lman_volume else 0.0

    tallied_results = {
        "Total Volume of Area X (um^3)": total_area_x_volume,
        "Volume of Lesion Overlap with Area X (um^3)": total_overlap_volume_x,
        "Percent of Area X Lesioned (%)": percent_area_x_lesioned,
        "Total Volume of LMAN (um^3)": total_lman_volume,
        "Volume of Lesion Overlap with LMAN (um^3)": total_overlap_volume_lman,
        "Percent of LMAN Lesioned (%)": percent_lman_lesioned,
        "Total Volume of Lesion Outside Striatum (Anterior) (um^3)": total_lesion_outside_striatum,
    }

    print(f"Reading file: {results_csv_path}")
    for k, v in tallied_results.items():
        print(f"{k}: {v}")
    return tallied_results


def quantify_section_areas_and_volumes(cfg: SeriesQuantConfig) -> QuantifyOutput:
    """
    Quantify areas (px², μm²) and volumes (μm³) for Lesion, Area X, LMAN across all rows.
    If the lesion intersects the striatum line, split the lesion polygon with that line,
    show a labeled plot with 'Lesion Part 1' and 'Lesion Part 2', and prompt the user
    to enter **1 or 2** for which part lies within the striatum.

    Always writes:
      - per-row CSV:  '<stem>_areas_and_volumes.csv'
      - summary JSON: '<stem>_final_volumes.json' (if cfg.write_summary_json)

    Files are saved in the chosen output_dir (asked interactively if not provided).
    """
    csv_path = Path(cfg.csv_path)
    df = pd.read_csv(csv_path).copy()

    # Build Section Number from 'Image Name' if present; otherwise keep None
    if "Image Name" in df.columns:
        names = df["Image Name"].astype(str)
    else:
        names = pd.Series([None] * len(df), index=df.index)

    df["Section Number"] = names.apply(lambda s: _get_section_number(s or "", cfg.section_regex))

    # Only sort if any section numbers are present
    if df["Section Number"].notna().any():
        df = df.sort_values(["Section Number"]).reset_index(drop=True)
    else:
        df = df.reset_index(drop=True)

    rows: List[Dict[str, Union[str, int, float, None]]] = []

    um_per_px = float(cfg.microns_per_pixel)
    area_scale = um_per_px ** 2  # px² -> μm²

    for idx, row in df.iterrows():
        # Parse coords
        str_x = _parse_coords_cell(row.get("Striatum_x_coords"))
        str_y = _parse_coords_cell(row.get("Striatum_y_coords"))
        les_x = _parse_coords_cell(row.get("Lesion_x_coords"))
        les_y = _parse_coords_cell(row.get("Lesion_y_coords"))
        axx_x = _parse_coords_cell(row.get("Area_X_x_coords"))
        axx_y = _parse_coords_cell(row.get("Area_X_y_coords"))
        lmn_x = _parse_coords_cell(row.get("LMAN_x_coords"))
        lmn_y = _parse_coords_cell(row.get("LMAN_y_coords"))

        lesion_poly = _coords_to_polygon(les_x, les_y)
        area_x_poly = _coords_to_polygon(axx_x, axx_y)
        lman_poly   = _coords_to_polygon(lmn_x, lmn_y)
        str_line    = _coords_to_line(str_x, str_y)

        # Areas in px²
        lesion_area_px2   = lesion_poly.area if lesion_poly else 0.0
        area_x_area_px2   = area_x_poly.area if area_x_poly else 0.0
        lman_area_px2     = lman_poly.area if lman_poly else 0.0
        overlap_x_px2     = lesion_poly.intersection(area_x_poly).area if (lesion_poly and area_x_poly and lesion_poly.intersects(area_x_poly)) else 0.0
        overlap_lman_px2  = lesion_poly.intersection(lman_poly).area if (lesion_poly and lman_poly and lesion_poly.intersects(lman_poly)) else 0.0

        # Convert to μm²
        lesion_area_um2  = lesion_area_px2 * area_scale
        area_x_area_um2  = area_x_area_px2 * area_scale
        lman_area_um2    = lman_area_px2 * area_scale
        overlap_x_um2    = overlap_x_px2 * area_scale
        overlap_lman_um2 = overlap_lman_px2 * area_scale

        # Effective thickness by section gaps (if any)
        sect_num = row.get("Section Number")
        if idx == 0 or pd.isna(sect_num):
            eff_thick_um = cfg.section_thickness_um
        else:
            prev_sect = df.loc[idx - 1, "Section Number"]
            if pd.isna(prev_sect):
                eff_thick_um = cfg.section_thickness_um
            else:
                gap = max(1, int(sect_num) - int(prev_sect))
                eff_thick_um = gap * cfg.section_thickness_um

        # Volumes (μm³)
        lesion_vol_um3  = lesion_area_um2  * eff_thick_um
        area_x_vol_um3  = area_x_area_um2  * eff_thick_um
        lman_vol_um3    = lman_area_um2    * eff_thick_um
        overlap_x_um3   = overlap_x_um2    * eff_thick_um
        overlap_lman_um3= overlap_lman_um2 * eff_thick_um

        # Split lesion by striatum line (interactive 1/2 choice)
        inside_um2 = outside_um2 = inside_um3 = outside_um3 = 0.0
        if lesion_poly and str_line and lesion_poly.intersects(str_line):
            parts = shapely_split(lesion_poly, str_line)
            if not parts.is_empty and len(parts.geoms) >= 2:
                # Take two largest pieces (robust vs tiny slivers)
                parts_sorted = sorted(parts.geoms, key=lambda p: p.area, reverse=True)
                p1, p2 = parts_sorted[0], parts_sorted[1]

                a1_um2 = p1.area * area_scale
                a2_um2 = p2.area * area_scale
                v1_um3 = a1_um2 * eff_thick_um
                v2_um3 = a2_um2 * eff_thick_um

                # Plot to help choose
                if cfg.plot_split_when_prompt:
                    plt.figure(figsize=(7, 7))
                    xs, ys = str_line.xy
                    plt.plot(xs, ys, color="tab:blue", lw=1.8, label="Striatum (Line)")
                    lx, ly = lesion_poly.exterior.xy
                    plt.fill(lx, ly, color="tab:red", alpha=0.25, label="Lesion Area")
                    for i, part in enumerate([p1, p2], start=1):
                        px, py = part.exterior.xy
                        plt.fill(px, py, alpha=0.5, label=f"Lesion Part {i}")
                        cx, cy = part.centroid.coords[0]
                        plt.text(cx, cy, str(i), ha="center", va="center", fontsize=12)
                    plt.xlabel("X Coordinates"); plt.ylabel("Y Coordinates")
                    plt.title(f"Lesion Area Split by Striatum for {row.get('Image Name', f'row{idx}')}")
                    plt.legend(); plt.gca().invert_yaxis()
                    plt.gca().set_aspect("equal", adjustable="box")
                    plt.tight_layout(); plt.show()

                # Console prompt: strictly 1 or 2
                ans = input(f"[{row.get('Image Name','row'+str(idx))}] Enter the part number (1 or 2) that is within the striatum: ").strip()
                if ans == "1":
                    inside_um2, outside_um2 = a1_um2, a2_um2
                    inside_um3, outside_um3 = v1_um3, v2_um3
                elif ans == "2":
                    inside_um2, outside_um2 = a2_um2, a1_um2
                    inside_um3, outside_um3 = v2_um3, v1_um3
                else:
                    print("Invalid input (expected '1' or '2'). Inside/outside set to 0 for this row.")

        rows.append({
            "Image Name": row.get("Image Name"),
            "Section Number": sect_num,
            "Effective Thickness (um)": eff_thick_um,

            "Lesion Area (px^2)": lesion_area_px2,
            "Area X Area (px^2)": area_x_area_px2,
            "LMAN Area (px^2)":   lman_area_px2,
            "Overlap Area with Area X (px^2)": overlap_x_px2,
            "Overlap Area with LMAN (px^2)":   overlap_lman_px2,

            "Lesion Area (um^2)": lesion_area_um2,
            "Area X Area (um^2)": area_x_area_um2,
            "LMAN Area (um^2)":   lman_area_um2,
            "Overlap Area with Area X (um^2)": overlap_x_um2,
            "Overlap Area with LMAN (um^2)":   overlap_lman_um2,

            "Lesion Volume (um^3)": lesion_vol_um3,
            "Area X Volume (um^3)": area_x_vol_um3,
            "LMAN Volume (um^3)":   lman_vol_um3,
            "Overlap Volume with Area X (um^3)": overlap_x_um3,
            "Overlap Volume with LMAN (um^3)":   overlap_lman_um3,

            "Lesion Area inside Striatum (um^2)": inside_um2,
            "Lesion Volume inside Striatum (um^3)": inside_um3,
            "Lesion Area outside Striatum (um^2)": outside_um2,
            "Lesion Volume outside Striatum (um^3)": outside_um3,
        })

    out_df = pd.DataFrame(rows)

    # Resolve and create output_dir (prompt if needed)
    out_dir = _choose_output_dir(cfg)
    base = Path(cfg.csv_path).with_suffix("").name

    # Always write per-row CSV in output_dir
    per_row_csv = out_dir / f"{base}_areas_and_volumes.csv"
    out_df.to_csv(per_row_csv, index=False)
    print("Saved:", per_row_csv)

    # Compute tallies and print
    tallies = _tally_totals(out_df, per_row_csv)

    # Optionally write the summary JSON in output_dir
    summary_json: Optional[Path] = None
    if cfg.write_summary_json:
        summary_json = out_dir / f"{base}_final_volumes.json"
        with open(summary_json, "w") as f:
            json.dump(tallies, f, indent=4)
        print(f"Analysis complete. Results saved to '{summary_json}'")

    return QuantifyOutput(
        dataframe=out_df,
        per_row_csv=per_row_csv,
        tallies=tallies,
        summary_json=summary_json,
    )


"""
from quantify_section_areas_and_volumes import SeriesQuantConfig, quantify_section_areas_and_volumes

csv_path = "/Users/mirandahulsey-vincent/Documents/allPythonCode/Histology_analysis/inputs/USA5510_Lower_Hemisphere_Lesion_Area_converted_for_validation.csv"

um_per_px = 568.0 / 500.0  # if 500 px == 568 μm

cfg = SeriesQuantConfig(
    csv_path=csv_path,
    microns_per_pixel=um_per_px,
    section_thickness_um=75.0,
    section_regex=r"sect(\d+)",
    plot_split_when_prompt=True,
    output_dir= "/Users/mirandahulsey-vincent/Documents/allPythonCode/Histology_analysis/inputs"
)

res = quantify_section_areas_and_volumes(cfg)


"""