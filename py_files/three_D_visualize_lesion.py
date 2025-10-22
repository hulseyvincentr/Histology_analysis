# -*- coding: utf-8 -*-
# three_D_visualize_lesion.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union

import ast
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Shapely is used to validate/normalize polygons and provide closed rings
try:
    from shapely.geometry import Polygon
except Exception as e:
    raise ImportError("This module requires 'shapely'. Install with: pip install shapely") from e


@dataclass
class Viz3DResult:
    fig: plt.Figure
    ax: plt.Axes
    saved_path: Path


# ----------------------------- Helpers -------------------------------------- #
_NONE_TOKENS = {
    "none", "none present", "none_present", "no shape", "not present",
    "nan", "", "0", "0.0"
}

def _parse_coords_cell(cell: Union[str, float, int, None]) -> Optional[List[float]]:
    """
    Accepts tuple-style '(1, 2, 3)', list '[1,2,3]', raw '1, 2, 3,',
    numbers (0/0.0 => None), or 'none present' tokens.
    Returns list[float] or None.
    """
    if cell is None:
        return None
    if isinstance(cell, (int, float)):
        return None if float(cell) == 0.0 else [float(cell)]
    s = str(cell).strip()
    if s.lower() in _NONE_TOKENS:
        return None
    # Try literal first
    try:
        obj = ast.literal_eval(s)
        if isinstance(obj, (list, tuple)):
            out: List[float] = []
            for v in obj:
                try:
                    out.append(float(v))
                except Exception:
                    pass
            return out or None
    except Exception:
        pass
    # Fallback: strip parens/brackets & split by commas
    s = s.strip("()[]")
    parts = [p.strip() for p in s.split(",") if p.strip()]
    out: List[float] = []
    for p in parts:
        try:
            out.append(float(p))
        except Exception:
            pass
    return out or None


def _plot_prism(
    ax: plt.Axes,
    x: List[float],
    y: List[float],
    z_base: float,
    thickness_um: float,
    *,
    color: str,
    alpha: float,
    label: Optional[str] = None,
):
    """Draw a vertical prism for a polygon (top, bottom, and side faces)."""
    if len(x) < 3 or len(y) < 3:
        return
    poly = Polygon(list(zip(x, y)))
    if poly.is_empty or not poly.is_valid or poly.area <= 0:
        return
    rx, ry = poly.exterior.xy  # closed ring; last == first
    rx = list(rx)
    ry = list(ry)
    z_top = z_base + thickness_um

    # Bottom & top faces
    ax.add_collection3d(
        Poly3DCollection([list(zip(rx, ry, [z_base] * len(rx)))], color=color, alpha=alpha, label=label)
    )
    ax.add_collection3d(
        Poly3DCollection([list(zip(rx, ry, [z_top] * len(rx)))], color=color, alpha=alpha)
    )

    # Side faces
    for i in range(len(rx) - 1):
        verts = [
            (rx[i],     ry[i],     z_base),
            (rx[i + 1], ry[i + 1], z_base),
            (rx[i + 1], ry[i + 1], z_top),
            (rx[i],     ry[i],     z_top),
        ]
        ax.add_collection3d(Poly3DCollection([verts], color=color, alpha=alpha))


def _resolve_results_csv(results_path: Path, csv_path: Path) -> Path:
    """
    Accept a CSV path OR a directory and return the matching *_areas_and_volumes.csv.
    Preference:
      1) <results_dir>/<coords_stem>_areas_and_volumes.csv
      2) If not found, single *_areas_and_volumes.csv in the directory.
    """
    if results_path.is_file():
        return results_path
    if results_path.is_dir():
        expected = results_path / f"{csv_path.stem}_areas_and_volumes.csv"
        if expected.is_file():
            return expected
        candidates = sorted(results_path.glob("*_areas_and_volumes.csv"))
        if not candidates:
            raise FileNotFoundError(
                f"No '*_areas_and_volumes.csv' found in directory: {results_path}"
            )
        if len(candidates) > 1:
            names = "\n  - ".join(c.name for c in candidates)
            raise FileExistsError(
                "Multiple results files found; please pass the exact CSV path:\n"
                f"  - {names}"
            )
        return candidates[0]
    raise FileNotFoundError(f"results_path not found: {results_path}")


# ----------------------------- Public API ----------------------------------- #
def visualize_lesion_3d(
    csv_path: Union[str, Path],
    results_path: Union[str, Path],
    *,
    microns_per_pixel: Optional[float] = None,   # X/Y scaling (px -> μm) if provided
    dpi: int = 200,
    show: bool = True,
    alpha: float = 0.5,
    colors: Optional[Dict[str, str]] = None,     # keys: "Area_X", "LMAN", "Lesion"
    title: Optional[str] = None,
    invert_y: bool = False,                      # True to mimic image coordinates (Y down)
) -> Viz3DResult:
    """
    Build a stacked 3D visualization of Area X, LMAN, and Lesion prisms across sections,
    then save the PNG **into the results_path location** (same folder as the results CSV).

    Parameters
    ----------
    csv_path : str | Path
        Coordinate CSV with columns like:
        ['Image Name', 'Area_X_x_coords', 'Area_X_y_coords',
         'LMAN_x_coords', 'LMAN_y_coords', 'Lesion_x_coords', 'Lesion_y_coords', ...]

    results_path : str | Path
        Either:
          - the per-row results CSV (containing 'Effective Thickness (um)'), or
          - a directory containing '<coords_stem>_areas_and_volumes.csv' (preferred),
            or only one '*_areas_and_volumes.csv' file.

    microns_per_pixel : float, optional
        If provided (e.g., 568/500 when 500 px == 568 μm), X/Y will be scaled to μm.
        Z uses 'Effective Thickness (um)' directly.

    Returns
    -------
    Viz3DResult(fig, ax, saved_path)
    """
    csv_path = Path(csv_path)
    results_path = _resolve_results_csv(Path(results_path), csv_path)

    # Load data
    df = pd.read_csv(csv_path).copy()
    results_df = pd.read_csv(results_path).copy()
    if "Effective Thickness (um)" not in results_df.columns:
        raise ValueError("results_df is missing column 'Effective Thickness (um)'")
    print("Using results CSV:", results_path)

    # Align rows: prefer 'Image Name', else fallback to index alignment
    if "Image Name" in df.columns and "Image Name" in results_df.columns:
        merged = pd.merge(
            df,
            results_df[["Image Name", "Effective Thickness (um)"]],
            on="Image Name",
            how="left",
            validate="many_to_one",
        )
        if merged["Effective Thickness (um)"].isna().any():
            missing = merged.loc[merged["Effective Thickness (um)"].isna(), "Image Name"].tolist()
            raise ValueError(
                "Could not find 'Effective Thickness (um)' for some Image Names. "
                f"Missing examples: {missing[:5]}"
            )
        df = merged
    else:
        if len(df) != len(results_df):
            raise ValueError(
                "Cannot align by Image Name and lengths differ. "
                "Either include 'Image Name' in both files or ensure equal row counts."
            )
        df["Effective Thickness (um)"] = results_df["Effective Thickness (um)"].values

    # Colors
    default_colors = {"Area_X": "tab:green", "LMAN": "tab:purple", "Lesion": "tab:red"}
    colors = {**default_colors, **(colors or {})}

    # Setup plot
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Iterate rows and stack prisms
    z_base = 0.0
    xs_all: List[float] = []
    ys_all: List[float] = []
    zs_all: List[float] = []

    def _scaled(arr: Optional[List[float]]) -> Optional[List[float]]:
        if arr is None:
            return None
        return [a * microns_per_pixel for a in arr] if microns_per_pixel else arr

    for idx, row in df.iterrows():
        eff_th = float(row["Effective Thickness (um)"])

        axx_x = _scaled(_parse_coords_cell(row.get("Area_X_x_coords")))
        axx_y = _scaled(_parse_coords_cell(row.get("Area_X_y_coords")))
        lmn_x = _scaled(_parse_coords_cell(row.get("LMAN_x_coords")))
        lmn_y = _scaled(_parse_coords_cell(row.get("LMAN_y_coords")))
        les_x = _scaled(_parse_coords_cell(row.get("Lesion_x_coords")))
        les_y = _scaled(_parse_coords_cell(row.get("Lesion_y_coords")))

        label_area_x = "Area X" if idx == 0 else None
        label_lman   = "LMAN"   if idx == 0 else None
        label_lesion = "Lesion" if idx == 0 else None

        if axx_x and axx_y:
            _plot_prism(ax, axx_x, axx_y, z_base, eff_th, color=colors["Area_X"], alpha=alpha, label=label_area_x)
            xs_all.extend(axx_x); ys_all.extend(axx_y)
        if lmn_x and lmn_y:
            _plot_prism(ax, lmn_x, lmn_y, z_base, eff_th, color=colors["LMAN"], alpha=alpha, label=label_lman)
            xs_all.extend(lmn_x); ys_all.extend(lmn_y)
        if les_x and les_y:
            _plot_prism(ax, les_x, les_y, z_base, eff_th, color=colors["Lesion"], alpha=alpha, label=label_lesion)
            xs_all.extend(les_x); ys_all.extend(les_y)

        zs_all.extend([z_base, z_base + eff_th])
        z_base += eff_th

    # Labels & title
    unit_xy = "μm" if microns_per_pixel else "px"
    ax.set_xlabel(f"X ({unit_xy})")
    ax.set_ylabel(f"Y ({unit_xy})")
    ax.set_zlabel("Depth (μm)")
    if title is None:
        title = f"3D Visualization of Area X, LMAN, and Lesion — {csv_path.stem}"
    ax.set_title(title)

    if invert_y:
        ax.invert_yaxis()

    # Limits with padding
    if xs_all and ys_all and zs_all:
        pad_x = max(10.0, 0.05 * (max(xs_all) - min(xs_all)))
        pad_y = max(10.0, 0.05 * (max(ys_all) - min(ys_all)))
        ax.set_xlim(min(xs_all) - pad_x, max(xs_all) + pad_x)
        ax.set_ylim(min(ys_all) - pad_y, max(ys_all) + pad_y)
        ax.set_zlim(0, max(zs_all) * 1.05)

    # Deduplicate legend entries
    handles, labels = ax.get_legend_handles_labels()
    if labels:
        seen = {}
        dedup_h, dedup_l = [], []
        for h, l in zip(handles, labels):
            if l and l not in seen:
                seen[l] = True
                dedup_h.append(h); dedup_l.append(l)
        if dedup_l:
            ax.legend(dedup_h, dedup_l, loc="best")

    fig.tight_layout()

    # Save the figure into the SAME DIRECTORY AS results CSV
    save_dir = results_path.parent
    save_dir.mkdir(parents=True, exist_ok=True)
    out_path = save_dir / f"{csv_path.stem}_3d_regions.png"
    fig.savefig(out_path, dpi=dpi)

    if show:
        plt.show()
    else:
        plt.close(fig)

    print("Saved figure:", out_path)
    return Viz3DResult(fig=fig, ax=ax, saved_path=out_path)


"""

from three_D_visualize_lesion import visualize_lesion_3d

csv_path     = "/Users/mirandahulsey-vincent/Documents/allPythonCode/Histology_analysis/inputs/USA5510_Lower_Hemisphere_Lesion_Area_converted_for_validation.csv"
results_path = "/Users/mirandahulsey-vincent/Documents/allPythonCode/Histology_analysis/inputs"

um_per_px = 568.0 / 500.0  # if 500 px == 568 μm

res = visualize_lesion_3d(
    csv_path=csv_path,
    results_path=results_path,         # file OR directory is OK
    microns_per_pixel=um_per_px,       # omit if X/Y already in μm
    show=True,
    alpha=0.5,
    invert_y=False,
)
print("Saved figure:", res.saved_path)



"""