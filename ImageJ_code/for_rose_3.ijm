
//Enter the path to images here
directory = "/Volumes/my_own_ssd/microscope_images/USA5510/lower_hemisphere";
//Enter the path to a folder to save ROI files
ROI_folder_path  = "/Volumes/my_own_ssd/microscope_images/USA5510/lower_hemisphere/ROIs"

//Create table, set settings for dialog box, make list of files in image folder
Table.create("Lesion Area");
dialog_labels = newArray("Striatum_", "Lesion_", "Area_X_", "LMAN_");
dialog_defaults = newArray(1,1,1,1);
filelist = getFileList(directory)


function process_region(region_name, present_or_absent, tool_for_annotation, rowIndex) { 
	//If the user marked the region as absent, record that, and move to next region
		if (present_or_absent==0) {
			xpoints = "none present";
			ypoints = "none present";
		}
	//Otherwise, have them annotate the region, and save the ROI, and save the coorinates into the table
		if (present_or_absent==1) {
			setTool(tool_for_annotation);
			waitForUser("Annotate the "+region_name+" then press ok");
			roiManager("add");
			roiManager("Save", ROI_folder_path+base_image+"_"+region_name+".roi");
			getSelectionCoordinates(xpoints, ypoints);
			xpoints = String.join(xpoints);
			ypoints = String.join(ypoints);
			
		}
		Table.set(region_name+"_x_coords", rowIndex, xpoints);
		Table.set(region_name+"_y_coords", rowIndex, ypoints);	
		run("Select None");
}

//Loop through each file ending in .jpg
for (i = 0; i < lengthOf(filelist); i++) {
    if (endsWith(filelist[i], ".jpg")) { 
    	roiManager("reset");
        open(directory + File.separator + filelist[i]);
		
		base_image = getTitle();
		base_image = substring(base_image, 1, base_image.length);
		Table.set("Image Name", i, base_image);
		
		//Have the user mark what regions are present vs absent
		Dialog.create("Please uncheck any that are missing");
		Dialog.addCheckboxGroup(4, 1, dialog_labels, dialog_defaults);
		Dialog.show();
		//Record their inputs
		STR_present = Dialog.getCheckbox();
		LES_present = Dialog.getCheckbox();
		ARX_present = Dialog.getCheckbox();
		LMAN_present = Dialog.getCheckbox();

		//Process each region
		process_region("Striatum", STR_present, "freeline", i);
		process_region("Lesion", LES_present, "polygon", i);
		process_region("Area_X", ARX_present, "polygon", i);
		process_region("LMAN", LMAN_present, "polygon", i);
		
      	close("*");  
    } 
}








