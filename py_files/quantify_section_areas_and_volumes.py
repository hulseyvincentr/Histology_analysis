# -*- coding: utf-8 -*-
# quantify_section_areas_and_volumes.py
from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

import re
import math
import pandas as pd
import matplotlib.pyplot as plt

# Shapely (install if missing): pip install shapely
try:
    from shapely.geometry import Polygon, LineString
    from shapely.ops import split as shapely_split
except Exception as e:
    raise ImportError(
        "This module requires 'shapely'. Install with: pip install shapely"
    ) from e


# ──────────────────────────────────────────────────────────────────────────────
# Parsing
# ──────────────────────────────────────────────────────────────────────────────
_NONE_TOKENS = {
    "none", "none present", "none_present", "no shape", "not present",
    "na", "n/a", "nan", "", "0", "0.0"
}

def _parse_coords_cell(cell: Union[str, float, int, None]) -> Optional[List[float]]:
    """Robustly parse a coord cell into a list[float] or None."""
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


def _get_section_number(image_name: str, *, pattern: str = r"sect(\d+)") -> Optional[int]:
    m = re.search(pattern, str(image_name))
    return int(m.group(1)) if m else None


# ──────────────────────────────────────────────────────────────────────────────
# Core API
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class SeriesQuantConfig:
    csv_path: Union[str, Path]
    microns_per_pixel: float            # e.g., if 500 px = 568 μm, set 568/500
    section_thickness_um: float = 75.0  # thickness of one section in μm
    section_regex: str = r"sect(\d+)"   # how to read section number from Image Name
    choose_inside: str = "prompt"       # 'prompt' | 'overlap_area_x' | 'skip'
    plot_split_when_prompt: bool = True # show the split plot when prompting
    out_csv_path: Optional[Union[str, Path]] = None  # write a CSV if given


def quantify_section_areas_and_volumes(cfg: SeriesQuantConfig) -> pd.DataFrame:
    """
    Quantify areas (px², μm²) and volumes (μm³) across all rows in a CSV,
    compute overlaps (Lesion∩Area X, Lesion∩LMAN), and optionally split the
    Lesion by the Striatum line to report inside/outside lesion area/volume.

    Returns a DataFrame and optionally writes CSV to cfg.out_csv_path.
    """
    csv_path = Path(cfg.csv_path)
    df = pd.read_csv(csv_path).copy()

    # Sort by section number (handles gaps by multiplying thickness by the gap)
    df["Section Number"] = df.get("Image Name", "").apply(
        lambda s: _get_section_number(s, pattern=cfg.section_regex)
    )
    df = df.sort_values(["Section Number"]).reset_index(drop=True)

    rows: List[Dict[str, Union[str, int, float, None]]] = []

    # Precompute scaling factors
    μm_per_px = float(cfg.microns_per_pixel)
    area_scale = μm_per_px ** 2  # px² -> μm²

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

        # Build geometries
        lesion_poly = _coords_to_polygon(les_x, les_y)
        area_x_poly = _coords_to_polygon(axx_x, axx_y)
        lman_poly   = _coords_to_polygon(lmn_x, lmn_y)
        str_line    = _coords_to_line(str_x, str_y)

        # Areas (px²)
        lesion_area_px2   = lesion_poly.area if lesion_poly else 0.0
        area_x_area_px2   = area_x_poly.area if area_x_poly else 0.0
        lman_area_px2     = lman_poly.area if lman_poly else 0.0
        overlap_x_px2     = lesion_poly.intersection(area_x_poly).area if (lesion_poly and area_x_poly and lesion_poly.intersects(area_x_poly)) else 0.0
        overlap_lman_px2  = lesion_poly.intersection(lman_poly).area if (lesion_poly and lman_poly and lesion_poly.intersects(lman_poly)) else 0.0

        # Areas (μm²)
        lesion_area_um2  = lesion_area_px2 * area_scale
        area_x_area_um2  = area_x_area_px2 * area_scale
        lman_area_um2    = lman_area_px2 * area_scale
        overlap_x_um2    = overlap_x_px2 * area_scale
        overlap_lman_um2 = overlap_lman_px2 * area_scale

        # Effective thickness by section gap
        sect_num = row.get("Section Number")
        if idx == 0 or pd.isna(sect_num):
            eff_thick_um = cfg.section_thickness_um
        else:
            prev_sect = df.loc[idx - 1, "Section Number"]
            if pd.isna(prev_sect):
                eff_thick_um = cfg.section_thickness_um
            else:
                gap = max(1, int(sect_num) - int(prev_sect))  # at least 1
                eff_thick_um = gap * cfg.section_thickness_um

        # Volumes (μm³) = area (μm²) × thickness (μm)
        lesion_vol_um3  = lesion_area_um2  * eff_thick_um
        area_x_vol_um3  = area_x_area_um2  * eff_thick_um
        lman_vol_um3    = lman_area_um2    * eff_thick_um
        overlap_x_um3   = overlap_x_um2    * eff_thick_um
        overlap_lman_um3= overlap_lman_um2 * eff_thick_um

        # Split lesion by striatum
        inside_um2 = outside_um2 = inside_um3 = outside_um3 = 0.0
        if lesion_poly and str_line and lesion_poly.intersects(str_line):
            parts = shapely_split(lesion_poly, str_line)
            if not parts.is_empty and len(parts.geoms) >= 2:
                p1, p2 = parts.geoms[:2]
                a1_um2 = p1.area * area_scale
                a2_um2 = p2.area * area_scale
                v1_um3 = a1_um2 * eff_thick_um
                v2_um3 = a2_um2 * eff_thick_um

                choice = None
                if cfg.choose_inside == "overlap_area_x" and area_x_poly:
                    # Choose the part with larger overlap with Area X
                    o1 = p1.intersection(area_x_poly).area
                    o2 = p2.intersection(area_x_poly).area
                    choice = 1 if o1 >= o2 else 2
                elif cfg.choose_inside == "prompt":
                    # Quick plot to let the user pick
                    if cfg.plot_split_when_prompt:
                        plt.figure(figsize=(7, 7))
                        # Striatum line
                        if str_line:
                            xs, ys = str_line.xy
                            plt.plot(xs, ys, color="tab:blue", lw=1.6, label="Striatum")
                        # Original lesion
                        lx, ly = lesion_poly.exterior.xy
                        plt.fill(lx, ly, color="tab:red", alpha=0.25, label="Lesion")
                        # Parts
                        for i, part in enumerate([p1, p2], start=1):
                            px, py = part.exterior.xy
                            plt.fill(px, py, alpha=0.5, label=f"Part {i}")
                        plt.gca().invert_yaxis()
                        plt.gca().set_aspect("equal", adjustable="box")
                        plt.legend()
                        plt.title(f"Split lesion — {row.get('Image Name', f'row{idx}')}")
                        plt.xlabel("X (px)"); plt.ylabel("Y (px)")
                        plt.tight_layout(); plt.show()

                    # Prompt
                    _ans = input(f"[{row.get('Image Name','row'+str(idx))}] Which part is inside Striatum? Enter 1 or 2: ").strip()
                    choice = 1 if _ans == "1" else 2 if _ans == "2" else None

                # assign inside/outside (default: larger area as 'inside' if no choice)
                if choice == 1:
                    inside_um2, outside_um2 = a1_um2, a2_um2
                    inside_um3, outside_um3 = v1_um3, v2_um3
                elif choice == 2:
                    inside_um2, outside_um2 = a2_um2, a1_um2
                    inside_um3, outside_um3 = v2_um3, v1_um3
                else:
                    # fallback: label the larger part as "inside"
                    if a1_um2 >= a2_um2:
                        inside_um2, outside_um2 = a1_um2, a2_um2
                        inside_um3, outside_um3 = v1_um3, v2_um3
                    else:
                        inside_um2, outside_um2 = a2_um2, a1_um2
                        inside_um3, outside_um3 = v2_um3, v1_um3

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

    # Optional write
    if cfg.out_csv_path:
        out_path = Path(cfg.out_csv_path)
        if out_path.is_dir():
            base = Path(cfg.csv_path).with_suffix("").name
            out_path = out_path / f"{base}_areas_and_volumes.csv"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(out_path, index=False)
        print("Saved:", out_path)

    return out_df


"""
from pathlib import Path
from quantify_section_areas_and_volumes import SeriesQuantConfig, quantify_section_areas_and_volumes

csv_path = "/Users/mirandahulsey-vincent/Documents/allPythonCode/Histology_analysis/inputs/USA5510_Lower_Hemisphere_Lesion_Area_converted_for_validation.csv"

# If 500 px == 568 μm, then microns_per_pixel = 568/500
μm_per_px = 568 / 500

cfg = SeriesQuantConfig(
    csv_path=csv_path,
    microns_per_pixel=μm_per_px,
    section_thickness_um=75.0,
    section_regex=r"sect(\d+)",          # adjust if your names differ
    choose_inside="prompt",               # 'prompt' | 'overlap_area_x' | 'skip'
    plot_split_when_prompt=True,
    out_csv_path=Path(csv_path).with_suffix("").parent,  # save alongside input
)

results_df = quantify_section_areas_and_volumes(cfg)
print(results_df.head(10))


"""