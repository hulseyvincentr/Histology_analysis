// ─────────────────────────────────────────────────────────────
// Paths (edit these):
// ─────────────────────────────────────────────────────────────
directory       = "/Volumes/my_own_ssd/microscope_images/USA5483/upper_hemisphere";
ROI_folder_path = "/Volumes/my_own_ssd/microscope_images/USA5483/upper_hemisphere_ROIs/";
table_save_path = ROI_folder_path + "Lesion_Area.csv";

// Regions and default tools (used only when we need to annotate)
region_names = newArray("Striatum", "Lesion", "Area_X", "LMAN");
region_tools = newArray("freeline", "polygon", "polygon", "polygon");

// Debug prints (set to 0 to silence)
DEBUG = 1;

// ─────────────────────────────────────────────────────────────
// Normalize paths so they always end with a separator
// ─────────────────────────────────────────────────────────────
sep = File.separator;
if (!endsWith(directory, sep))       directory       = directory + sep;
if (!endsWith(ROI_folder_path, sep)) ROI_folder_path = ROI_folder_path + sep;

// Make sure ROI folder exists
File.makeDirectory(ROI_folder_path);

// ─────────────────────────────────────────────────────────────
// Prep
// ─────────────────────────────────────────────────────────────
Table.create("Lesion Area");
filelist     = getFileList(directory);        // image filenames only
roi_filelist = getFileList(ROI_folder_path);  // ROI filenames only

// Membership check in roi_filelist
function roi_name_in_folder(roi_filename) {
    for (q=0; q<roi_filelist.length; q++) {
        if (roi_filelist[q] == roi_filename) return 1;
    }
    return 0;
}

// "<image>.jpg_<Region>.roi"
function roi_filename_for(base_image, region_name) {
    return base_image + "_" + region_name + ".roi";
}

// Count existing ROIs for an image; fill exists_flags[r] = 1/0
function roi_existence_by_listing(base_image, exists_flags) {
    count = 0;
    for (r=0; r<region_names.length; r++) {
        roi_fname = roi_filename_for(base_image, region_names[r]);
        if (roi_name_in_folder(roi_fname)==1) {
            exists_flags[r] = 1;
            count++;
        } else {
            exists_flags[r] = 0;
        }
    }
    return count;
}

// Record coordinates from current selection
function record_selection_coords(region_name, rowIndex) {
    getSelectionCoordinates(xpoints, ypoints);
    xpoints = String.join(xpoints);
    ypoints = String.join(ypoints);
    Table.set(region_name + "_x_coords", rowIndex, xpoints);
    Table.set(region_name + "_y_coords", rowIndex, ypoints);
    run("Select None");
}

// Load ROI and record coords (by filename membership)
function load_and_record_if_present(base_image, region_name, rowIndex) {
    roi_fname = roi_filename_for(base_image, region_name);
    if (roi_name_in_folder(roi_fname)==1) {
        fullpath = ROI_folder_path + roi_fname;  // safe: path normalized above
        roiManager("reset");
        roiManager("Open", fullpath);
        roiManager("Select", 0);
        record_selection_coords(region_name, rowIndex);
        return 1;
    } else {
        Table.set(region_name + "_x_coords", rowIndex, "none present");
        Table.set(region_name + "_y_coords", rowIndex, "none present");
        return 0;
    }
}

// Interactively annotate & save when present==1, else mark as absent
function annotate_or_absent(region_name, tool_for_annotation, present, rowIndex, base_image) {
    if (present==1) {
        setTool(tool_for_annotation);
        waitForUser("Annotate the " + region_name + " then press OK");
        roiManager("reset");
        roiManager("add");
        savepath = ROI_folder_path + roi_filename_for(base_image, region_name);
        roiManager("Save", savepath);

        // ── FIX: refresh the ROI file list from disk (avoids Array.concat on a string)
        roi_filelist = getFileList(ROI_folder_path);

        record_selection_coords(region_name, rowIndex);
    } else {
        Table.set(region_name + "_x_coords", rowIndex, "none present");
        Table.set(region_name + "_y_coords", rowIndex, "none present");
    }
}

// ─────────────────────────────────────────────────────────────
// Main loop
// ─────────────────────────────────────────────────────────────
row = 0;
for (i=0; i<lengthOf(filelist); i++) {
    if (!endsWith(filelist[i], ".jpg")) continue;

    base_image = filelist[i];  // e.g., "USA5510_..._sect115_lower.jpg"
    if (DEBUG && i<3) print("Image:", base_image);

    open(directory + base_image);
    Table.set("Image Name", row, base_image);

    exists_flags   = newArray(region_names.length);
    existing_count = roi_existence_by_listing(base_image, exists_flags);

    if (DEBUG && i<3) {
        for (r=0; r<region_names.length; r++) {
            print("  ROI present? ", region_names[r], " -> ", exists_flags[r]);
        }
        print("  Total existing for image: ", existing_count);
    }

    if (existing_count > 0) {
        // NON-INTERACTIVE: load what exists; mark others absent
        for (r=0; r<region_names.length; r++) {
            if (exists_flags[r]==1) {
                load_and_record_if_present(base_image, region_names[r], row);
            } else {
                Table.set(region_names[r] + "_x_coords", row, "none present");
                Table.set(region_names[r] + "_y_coords", row, "none present");
            }
        }
    } else {
        // INTERACTIVE: no ROIs at all → one presence dialog, then annotate/save
        Dialog.create("Annotate regions for this image (uncheck any that are absent)");
        Dialog.addCheckboxGroup(region_names.length, 1, region_names, newArray(1,1,1,1));
        Dialog.show();

        present_flags = newArray(region_names.length);
        for (r=0; r<region_names.length; r++) {
            present_flags[r] = Dialog.getCheckbox(); // 1 or 0
        }
        for (r=0; r<region_names.length; r++) {
            annotate_or_absent(region_names[r], region_tools[r], present_flags[r], row, base_image);
        }
    }

    close("*");
    row++;
}

// ─────────────────────────────────────────────────────────────
// Auto-save the table at the end
// ─────────────────────────────────────────────────────────────
if (isOpen("Lesion Area")) {
    selectWindow("Lesion Area");
    Table.save(table_save_path);
    print("Saved table to: " + table_save_path);
} else {
    print("Warning: 'Lesion Area' table not found; nothing saved.");
}
