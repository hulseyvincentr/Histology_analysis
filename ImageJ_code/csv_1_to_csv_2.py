from pathlib import Path
import pandas as pd, ast

in_path  = Path("/Volumes/my_own_ssd/microscope_images/USA5510/lower_hemisphere_ROIs/USA5510_Lower_Hemisphere_Lesion_Area.csv")
out_path = Path("/Volumes/my_own_ssd/microscope_images/USA5510/lower_hemisphere_ROIs/USA5510_Lower_Hemisphere_Lesion_Area_wide.csv")
keep_jpg_in_image_name = True

REGIONS = ["Striatum","Lesion","Area_X","LMAN"]

def _norm_region(value:str)->str:
    s = str(value).strip().replace(" ","_")
    if not s.endswith(".roi"): s += ".roi"
    return s

def _image_key(slide_name:str, keep_jpg:bool=True)->str:
    s = str(slide_name)
    if ".jpg" in s:
        base = s.split(".jpg")[0]
        return base + (".jpg" if keep_jpg else "")
    return s

def _coords_to_xy_strings(coords_str:str):
    if coords_str is None: return "none present","none present"
    txt = str(coords_str).strip()
    if txt.lower()=="none present" or txt=="": return "none present","none present"
    import ast
    try:
        parsed = ast.literal_eval(txt)
        xs, ys = zip(*parsed)
        def _s(v):
            try: return str(int(float(v)))
            except: return str(v)
        return ", ".join(_s(v) for v in xs), ", ".join(_s(v) for v in ys)
    except Exception:
        return "none present","none present"

def convert_long_to_wide(df_long, keep_jpg=True):
    required = {"slide_name","region","coordinates"}
    missing = required - set(df_long.columns)
    if missing:
        raise ValueError(f"Input CSV missing columns: {sorted(missing)}")
    df = df_long.copy()
    df["__image"]  = df["slide_name"].apply(lambda s: _image_key(s, keep_jpg))
    df["__region"] = df["region"].apply(_norm_region)
    rows = []
    for image, g in df.groupby("__image", sort=False):
        row = {"Image Name": image}
        for r in REGIONS:
            row[f"{r}_x_coords"] = "none present"
            row[f"{r}_y_coords"] = "none present"
        for _, rec in g.iterrows():
            base = str(rec["__region"]).replace(".roi","")
            if base in REGIONS:
                x_str, y_str = _coords_to_xy_strings(rec["coordinates"])
                row[f"{base}_x_coords"] = x_str
                row[f"{base}_y_coords"] = y_str
        rows.append(row)
    ordered = ["Image Name"] + [c for r in REGIONS for c in (f"{r}_x_coords", f"{r}_y_coords")]
    return pd.DataFrame(rows, columns=ordered)

df_long = pd.read_csv(in_path)
df_wide = convert_long_to_wide(df_long, keep_jpg=keep_jpg_in_image_name)
print(df_wide.head(10))
print("Shape:", df_wide.shape)
df_wide.to_csv(out_path, index=False)
print("Wrote:", out_path)