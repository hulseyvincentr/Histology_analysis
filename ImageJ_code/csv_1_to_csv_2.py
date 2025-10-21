# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 21:04:18 2024

@author: heino
"""

import pandas as pd
import ast

rose_dict = {"Image Name": [],
             "Striatum_x_coords":[],
             "Striatum_y_coords":[],
             "Lesion_x_coords":[],
             "Lesion_y_coords":[],
             "Area_X_x_coords":[],
             "Area_X_y_coords":[],
             "LMAN_x_coords":[],
             "LMAN_y_coords":[],
             }

df = pd.read_csv(r"""/Users/mirandahulsey-vincent/Documents/allPythonCode/Histology_analysis/inputs/sample_shapes_for_validation.csv""")

expected_regions = ["Area_X.roi", "Striatum.roi", "Lesion.roi", "LMAN.roi"]

unique_image_names = []
all_images_names = []
for i in range(len(df["slide_name"])):
    slide_name = df["slide_name"].iloc[i]
    image_name = slide_name.split(".jpg")[0]
    
    all_images_names.append(image_name)
    
df["image_names"] = all_images_names
    
unique_image_names = df["image_names"].unique()
    
for image_name in unique_image_names:
    img_filter = df["image_names"] == image_name
    df_img = df[img_filter]
    
    regions = list(df_img["region"])

    rose_dict["Image Name"].append(image_name)
    missing_items = list(set(expected_regions) - set(regions))

    for missing_item in missing_items:
        if missing_item == "Area_X.roi":
            rose_dict["Area_X_x_coords"].append("none present")
            rose_dict["Area_X_y_coords"].append("none present")
        if missing_item == "Striatum.roi":
            rose_dict["Striatum_x_coords"].append("none present")
            rose_dict["Striatum_y_coords"].append("none present")
        if missing_item == "Lesion.roi":
            rose_dict["Lesion_x_coords"].append("none present")
            rose_dict["Lesion_y_coords"].append("none present")
        if missing_item == "LMAN.roi":
            rose_dict["LMAN_x_coords"].append("none present")
            rose_dict["LMAN_y_coords"].append("none present")

    for i in range(len(df_img["slide_name"])):
        coords = df_img["coordinates"].iloc[i]
        region = df_img["region"].iloc[i]
        
        print(df_img)
        print(coords[0:10])
        print(region)
        result = ast.literal_eval(coords)

        
        list1, list2 = zip(*result)
        
        if region == "Area_X.roi":
            rose_dict["Area_X_x_coords"].append(list1)
            rose_dict["Area_X_y_coords"].append(list2)
        elif region == "Striatum.roi":
            rose_dict["Striatum_x_coords"].append(list1)
            rose_dict["Striatum_y_coords"].append(list2)
        elif region == "Lesion.roi":
            rose_dict["Lesion_x_coords"].append(list1)
            rose_dict["Lesion_y_coords"].append(list2)
        elif region == "LMAN.roi":
            rose_dict["LMAN_x_coords"].append(list1)
            rose_dict["LMAN_y_coords"].append(list2)
        

df_2 = pd.DataFrame(rose_dict)

df_2.to_csv(r"""/Users/mirandahulsey-vincent/Documents/allPythonCode/Histology_analysis/inputs/sample_shapes_for_validation.csv""")
    