# -*- coding: utf-8 -*-
# visualize_and_quantify_section_area.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import re
import math
import pandas as pd
import matplotlib.pyplot as plt


# ──────────────────────────────────────────────────────────────────────────────
# Public result container
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class SectionQuantResults:
    csv_path: Path
    row_index: int
    image_name: Optional[str]
    # Areas in pixel^2
    lesion_area_px2: float
    area_x_area_px2: float
    lman_area_px2: float
    # Optional metric for the striatum line (polyline length)
    striatum_length_px: Optional[float]
    # Scaled metrics (if scale_microns_per_pixel is provided)
    lesion_area_mm2: Optional[float] = None
    area_x_area_mm2: Optional[float] = None
    lman_area_mm2: Optional[float] = None
    striatum_length_mm: Optional[float] = None
    # Path to saved figure (if saved)
    figure_path: Optional[Path] = None


# ──────────────────────────────────────────────────────────────────────────────
# Robust parsing
# ──────────────────────────────────────────────────────────────────────────────
_NONE_TOKENS = {
    "none", "none present", "none_present", "no shape", "not present",
    "na", "n/a", "nan", "", "0", "0.0"
}

def _parse_coord_cell(cell: Union[str, float, int, None]) -> Optional[List[float]]:
    """
    Accepts:
      - tuple-like strings '(1, 2, 3)', CSV strings '1, 2, 3,', bracketed '[]'
      - raw comma lists '1,2,3'
      - sentinel values meaning 'no shape' (0.0, 'none present', etc.)
    Returns list[float] or None.
    """
    if cell is None:
        return None
    # Treat numeric 0 or 0.0 as 'none'
    if isinstance(cell, (int, float)) and float(cell) == 0.0:
        return None

    s = str(cell).strip()
    if s.lower() in _NONE_TOKENS:
        return None

    # normalize '(...)' / '[...]'
    s = s.strip("()[]")
    # split on commas, tolerate extra spaces & trailing comma
    parts = [p.strip() for p in s.split(",") if p.strip() != ""]
    out: List[float] = []
    for p in parts:
        # extract first number-like token per piece
        m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", p)
        if m:
            try:
                out.append(float(m.group(0)))
            except Exception:
                pass
    return out or None


# ──────────────────────────────────────────────────────────────────────────────
# Geometry helpers
# ──────────────────────────────────────────────────────────────────────────────
def _shoelace_area(x: Sequence[float], y: Sequence[float]) -> float:
    """Area via Shoelace formula (pixels^2). Auto-closes polygon if needed."""
    n = min(len(x), len(y))
    if n < 3:
        return 0.0
    # close polygon without modifying inputs
    xs = list(x[:n]) + [x[0]]
    ys = list(y[:n]) + [y[0]]
    s = 0.0
    for i in range(n):
        s += xs[i] * ys[i + 1] - ys[i] * xs[i + 1]
    return 0.5 * abs(s)


def _polyline_length(x: Sequence[float], y: Sequence[float]) -> float:
    """Length of a polyline (pixels)."""
    n = min(len(x), len(y))
    if n < 2:
        return 0.0
    total = 0.0
    for i in range(n - 1):
        dx = x[i + 1] - x[i]
        dy = y[i + 1] - y[i]
        total += math.hypot(dx, dy)
    return total


# ──────────────────────────────────────────────────────────────────────────────
# Core function
# ──────────────────────────────────────────────────────────────────────────────
def visualize_and_quantify_section_area(
    csv_path: Union[str, Path],
    *,
    row_index: Optional[int] = None,
    image_name: Optional[str] = None,
    scale_microns_per_pixel: Optional[float] = None,  # e.g., 0.5 μm/pixel
    save_path: Optional[Union[str, Path]] = None,
    show: bool = True,
    figsize: Tuple[float, float] = (8, 8),
    region_colors: Optional[Dict[str, str]] = None,
) -> SectionQuantResults:
    """
    Visualize and quantify traced regions from a single CSV row.

    Parameters
    ----------
    csv_path : str | Path
        CSV containing columns:
        ['Image Name', 'Striatum_x_coords', 'Striatum_y_coords', 'Lesion_x_coords',
         'Lesion_y_coords', 'Area_X_x_coords', 'Area_X_y_coords', 'LMAN_x_coords', 'LMAN_y_coords']
        Coordinate cells may be tuple-like strings '(…)', raw CSV '1,2,3', or 0.0/none tokens.

    row_index : int, optional
        Row index to visualize. Ignored if `image_name` is provided.
    image_name : str, optional
        Select row by exact 'Image Name' match (takes precedence over row_index).

    scale_microns_per_pixel : float, optional
        If provided, returns area in mm^2 and striatum length in mm.

    save_path : str | Path, optional
        If provided, saves the PNG here. If a directory is given, an auto
        filename '<image_name>_regions.png' is used.

    show : bool
        Whether to display the plot interactively.

    figsize : (w, h)
        Size for matplotlib figure.

    region_colors : dict, optional
        Override default colors: keys in {'Striatum','Lesion','Area_X','LMAN'}.

    Returns
    -------
    SectionQuantResults
        Areas (px^2 and optionally mm^2), striatum length, and saved figure path.
    """
    csv_path = Path(csv_path)
    df = pd.read_csv(csv_path)

    # choose row
    sel_idx: int
    if image_name is not None and "Image Name" in df.columns:
        matches = df.index[df["Image Name"] == image_name].tolist()
        if not matches:
            raise ValueError(f"image_name={image_name!r} not found in 'Image Name' column.")
        sel_idx = matches[0]
    else:
        if row_index is None:
            sel_idx = 0
        else:
            sel_idx = int(row_index)
        if sel_idx < 0 or sel_idx >= len(df):
            raise IndexError(f"row_index {sel_idx} out of range [0, {len(df)-1}]")

    row = df.iloc[sel_idx]
    img_name: Optional[str] = row["Image Name"] if "Image Name" in df.columns else None

    # parse coords
    str_x = _parse_coord_cell(row.get("Striatum_x_coords"))
    str_y = _parse_coord_cell(row.get("Striatum_y_coords"))
    les_x = _parse_coord_cell(row.get("Lesion_x_coords"))
    les_y = _parse_coord_cell(row.get("Lesion_y_coords"))
    axx_x = _parse_coord_cell(row.get("Area_X_x_coords"))
    axx_y = _parse_coord_cell(row.get("Area_X_y_coords"))
    lman_x = _parse_coord_cell(row.get("LMAN_x_coords"))
    lman_y = _parse_coord_cell(row.get("LMAN_y_coords"))

    # compute metrics (pixels)
    lesion_area_px2 = _shoelace_area(les_x or [], les_y or [])
    area_x_area_px2 = _shoelace_area(axx_x or [], axx_y or [])
    lman_area_px2   = _shoelace_area(lman_x or [], lman_y or [])
    str_len_px      = _polyline_length(str_x or [], str_y or []) if (str_x and str_y) else None

    # optional scaling to mm^2 / mm
    lesion_area_mm2 = area_x_area_mm2 = lman_area_mm2 = str_len_mm = None
    if scale_microns_per_pixel and scale_microns_per_pixel > 0:
        # area: (μm/pix)^2 -> μm^2 -> mm^2 (divide by 1e6)
        scale_area = (scale_microns_per_pixel ** 2) / 1_000_000.0
        lesion_area_mm2 = lesion_area_px2 * scale_area
        area_x_area_mm2 = area_x_area_px2 * scale_area
        lman_area_mm2   = lman_area_px2 * scale_area
        # length: μm/pix -> mm (divide by 1000)
        if str_len_px is not None:
            str_len_mm = (str_len_px * scale_microns_per_pixel) / 1000.0

    # plotting
    default_colors = {
        "Striatum": "tab:blue",
        "Lesion": "tab:red",
        "Area_X": "tab:green",
        "LMAN": "tab:purple",
    }
    colors = {**default_colors, **(region_colors or {})}

    fig, ax = plt.subplots(figsize=figsize)

    # Plot line and filled polygons
    if str_x and str_y:
        ax.plot(str_x, str_y, label=f"Striatum (line){'' if str_len_px is None else f' — len: {str_len_px:.1f}px'}",
                color=colors["Striatum"], linewidth=1.8)

    def _fill_if_any(x, y, label, color, area_px2: float, area_mm2: Optional[float]):
        if x and y and len(x) >= 3 and len(y) >= 3:
            pretty = f"{label} — area: {area_px2:.0f} px²"
            if area_mm2 is not None:
                pretty += f" ({area_mm2:.4f} mm²)"
            ax.fill(x, y, label=pretty, color=color, alpha=0.45, edgecolor=color, linewidth=1.0)

    _fill_if_any(les_x, les_y, "Lesion", colors["Lesion"], lesion_area_px2, lesion_area_mm2)
    _fill_if_any(axx_x, axx_y, "Area X", colors["Area_X"], area_x_area_px2, area_x_area_mm2)
    _fill_if_any(lman_x, lman_y, "LMAN", colors["LMAN"], lman_area_px2, lman_area_mm2)

    # axes cosmetics
    ax.set_xlabel("X (pixels)")
    ax.set_ylabel("Y (pixels)")
    title = f"Traced Regions{'' if img_name is None else f' — {img_name}'}"
    ax.set_title(title)
    ax.invert_yaxis()          # typical image coords
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    ax.legend(loc="best", frameon=True)

    # bounds
    all_x = []
    all_y = []
    for arr in (str_x, les_x, axx_x, lman_x):
        if arr: all_x.extend(arr)
    for arr in (str_y, les_y, axx_y, lman_y):
        if arr: all_y.extend(arr)
    if all_x and all_y:
        pad = max(10.0, 0.05 * max(max(all_x) - min(all_x), max(all_y) - min(all_y)))
        ax.set_xlim(min(all_x) - pad, max(all_x) + pad)
        ax.set_ylim(max(all_y) + pad, min(all_y) - pad)  # already inverted

    fig.tight_layout()

    # save/show
    fig_path: Optional[Path] = None
    if save_path is not None:
        save_path = Path(save_path)
        if save_path.is_dir():
            fname = f"{(img_name or f'row{sel_idx}').replace('/', '_')}_regions.png"
            fig_path = save_path / fname
        else:
            fig_path = save_path
        fig_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(fig_path, dpi=200)
    if show:
        plt.show()
    else:
        plt.close(fig)

    return SectionQuantResults(
        csv_path=csv_path,
        row_index=sel_idx,
        image_name=img_name,
        lesion_area_px2=lesion_area_px2,
        area_x_area_px2=area_x_area_px2,
        lman_area_px2=lman_area_px2,
        striatum_length_px=str_len_px,
        lesion_area_mm2=lesion_area_mm2,
        area_x_area_mm2=area_x_area_mm2,
        lman_area_mm2=lman_area_mm2,
        striatum_length_mm=str_len_mm,
        figure_path=fig_path,
    )


"""
###### Minimal: views on screen, no scaling
from visualize_and_quantify_section_area import visualize_and_quantify_section_area
from pathlib import Path

csv_path = "/Users/mirandahulsey-vincent/Documents/allPythonCode/Histology_analysis/inputs/USA5510_Lower_Hemisphere_Lesion_Area_converted_for_validation.csv"

res = visualize_and_quantify_section_area(
    csv_path,
    row_index=10,           # or image_name="USA5326_031524.02"
    show=True,             # display interactively
    save_path=None,        # or a directory/file path to save PNG
)

print(res)



############ With scale saves PNG
from visualize_and_quantify_section_area import visualize_and_quantify_section_area
from pathlib import Path

csv_path = "/Users/mirandahulsey-vincent/Documents/allPythonCode/Histology_analysis/inputs/Area_X_annotation_ROIs_and_csvs/USA5510_Lower_Hemisphere_Lesion_Area_converted_for_validation.csv"
outdir   = Path("/Users/mirandahulsey-vincent/Documents/allPythonCode/Histology_analysis/figures")
outdir.mkdir(parents=True, exist_ok=True)

res = visualize_and_quantify_section_area(
    csv_path,
    image_name=None,               # or provide exact 'Image Name'
    row_index=4,                   # used if image_name is None
    scale_microns_per_pixel=0.5,   # example pixel size; adjust to your imaging
    save_path=outdir,              # directory -> auto filename
    show=False,                    # just save
)

print("Lesion area:", res.lesion_area_px2, "px²", "/", res.lesion_area_mm2, "mm²")
print("Area X area:", res.area_x_area_px2, "px²", "/", res.area_x_area_mm2, "mm²")
print("LMAN area:",  res.lman_area_px2,   "px²", "/", res.lman_area_mm2,   "mm²")
print("Striatum length:", res.striatum_length_px, "px", "/", res.striatum_length_mm, "mm")
print("Saved figure:", res.figure_path)

"""