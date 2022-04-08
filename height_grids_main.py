"""
height_grids_main generate height grids from LAS point cloud data.
"""
import os
import shutil
import time

### IMPORT HEIGHT FUNCTIONS
from height_grids import get_merge_pipeline, \
    get_crop_pipeline, interpolation_grids
from utils import get_img_bbox

### PARAMETERS
from config import DIR_IMAGES_GEOTIFF, \
    LAS_INPUT_FOLDER, INTERPOLATION_FOLDER

######################### USER PARAMS #################################
### ENTER WHAT INTERPOLATION YOU WANT
no_interp= False
idw = False
nn = True
### TO GENERATE CROPPED POINT CLOUDS PER IMAGE WITHOUT INTERPOLATION YET, 
# ENTER TRUE FOR crop_las_files, DEFAULT IS TRUE FOR GENRATING INTERPOLATIONS
crop_las_files = False
### CLEAN-UP TO DELETE TEMP FILES CONTAINING PDAL .JSON PIPELINES
clean_up = True

# NECESSARY INPUT: A SINGLE MERGED POINT CLOUD OR PATH TO CREATE IT
las_merged = f"{LAS_INPUT_FOLDER}/merged.las"

########################################################################
# starting time
start = time.time()

for img_name in os.listdir(DIR_IMAGES_GEOTIFF)[:3]:
    if img_name[-4:] != ".tif":
        continue
    else:
        ### DEFINE ALL PATHS
        img_fp = DIR_IMAGES_GEOTIFF + "/" + img_name
        las_cropped = f"{LAS_INPUT_FOLDER}/buildings_cropped/{img_name[:-4]}.las"
 
        # GET IMAGE BBOX
        img_width, img_height, \
            [xmin, ymin, xmax, ymax], \
                [ulx_deg, uly_deg, lrx_deg, lry_deg] = get_img_bbox(img_fp)
        # bounding-box of the image needs to be in degrees !
        bbox_img = [(ulx_deg,lry_deg), (lrx_deg,uly_deg)]
        origin = [ulx_deg, lry_deg]

        ### IF BUILDING'S POINT CLOUD DOES NOT EXISTS
        if not os.path.exists(las_cropped) or crop_las_files == True:
            # WE MERGE ALL INPUT LAS IF THIS HAS NOT BEEN DONE YET
            if not os.path.exists(las_merged):
                print("No merged Point Cloud found... If you already had one, please check its name and directory.\
                    \nCreating a merged file for all LAS input.")
                in_srs, out_srs = 25832, 3857
                json_out_merge = f"{LAS_INPUT_FOLDER}/merge.json"
                # CREATE AND RUN PDAL MERGE-PIPELINE
                cmd_merge = get_merge_pipeline(LAS_INPUT_FOLDER, las_merged, in_srs, out_srs)
                with open (json_out_merge, "w") as outfile_merge:
                    outfile_merge.write(cmd_merge)
                os.system(f"pdal pipeline {json_out_merge}")
                os.remove(json_out_merge)
            cmd = get_crop_pipeline(xmin, xmax, ymin, ymax, las_merged, LAS_INPUT_FOLDER, img_name[:-4])
            if not os.path.exists(f"{LAS_INPUT_FOLDER}/buildings_cropped"):
                os.mkdir(f"{LAS_INPUT_FOLDER}/buildings_cropped")
                os.mkdir(f"{LAS_INPUT_FOLDER}/buildings_cropped/temp")
            json_out_bdg_crop = f"{LAS_INPUT_FOLDER}/buildings_cropped/temp/{img_name[:-4]}.json"
            # CREATE AND RUN PDAL CROPPING-PIPELINE
            with open (json_out_bdg_crop, "w") as outfile:
                outfile.write(cmd)
            os.system(f"pdal pipeline {json_out_bdg_crop}")
        # CREATE INTERPOLATION GRIDS
        interpolation_grids(img_fp, las_cropped, INTERPOLATION_FOLDER, \
            img_name, bbox_img, origin, \
            no_interp=no_interp, idw=idw, nn=nn)
    
# CLEAN UP TEMPORARY FOLDERS
if clean_up == True:
    temp_folder1 = f"{INTERPOLATION_FOLDER}/no_interpolation/temp"
    temp_folder2 = f"{INTERPOLATION_FOLDER}/idw_interpolation/temp"
    temp_folder3 = f"{LAS_INPUT_FOLDER}/buildings_cropped/temp"
    if os.path.exists(temp_folder3):
        shutil.rmtree(temp_folder3)
        print(f"Removed {temp_folder3}. Cropping the buildings point cloud is done.")
    if os.path.exists(temp_folder1):
        shutil.rmtree(temp_folder1)
        print(f"Removed {temp_folder1}. No-interpolation grid is finished.")
    if os.path.exists(temp_folder2):
        shutil.rmtree(temp_folder2)
        print(f"Removed {temp_folder2}. IDW-interpolation grid is finished.")

# end time
end = time.time()
print(f"Runtime of the program was {int((end - start)/60)} min, {int((end - start)%60)} sec.")