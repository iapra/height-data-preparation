import time
import numpy as np

from height_grids_main import height_grids_main
from relative_height_main import relative_height_main
from data_splitter import data_splitter
from utils import get_png, visualize_empty_height_data
from config import DIR_IMAGES_GEOTIFF, LAS_INPUT_FOLDER, DATA_SPLIT_FOLDER, \
    DIR_ABS_HEIGHT, DIR_REL_HEIGHT

# 1 # 
### UPDATE CONFIG FILE & DATA FOLDER

######################## RELATIVE HEIGHT ################################
# 2 #
### INDICATE THE LAS CLASSIFICATION YOU WANT
classification = 6
### TO GET NEW LAS FILE WITH RELATIVE HEIGHT, INDICATE TRUE AND THE DISTANCE TYPE
distance_type = str('vertical')    # 'vertical or 'min'
get_new_las = False

######################### INTERPOLATION #################################
# 3 #
### INDICATE A SINGLE MERGED POINT CLOUD OR A PATH TO CREATE IT
las_merged = f"{LAS_INPUT_FOLDER}/merged.las"
### ENTER WHAT INTERPOLATION YOU WANT WITH TRUE/FALSE
no_interp = False
idw = False
nn = False

######################### DATA-SPLIT ###################################
# 4 #
# INDICATE WHAT YOU NEED WITH TRUE/FALSE
create_split = False 
visualise_split = False
# IF visualise_split == True, give a directory:
DATA_SPLIT_FOLDER_TO_VISUALIZE = None

######################### CONVERT TIF TO PNG ###########################
# 5 #
### ENTER WHICH FOLDER(S) CONTAINING GEOTIFF IMAGES TO CONVERT & NODATA VALUE
FOLDER_LIST = []                    # leave list empty if no conversion needed
NODATA = -9999
dtype = np.uint8

#################### VISUALISE EMPTY HEIGHT DATA GRIDS #################
# 6 #
### ENTER WHICH FOLDER(S) CONTAINING PNG IMAGES TO ASSESS
HEIGHT_FOLDER_PATH = ''             # leave empty if no assessment needed
GEOJSON_EMPTY_TILES_PATH = ''       # leave empty if no assessment needed
NODATA_PNG = 0
########################################################################

if __name__ == "__main__":
    # starting time
    start = time.time()

    relative_height_main(distance_type, classification, get_new_las)

    height_grids_main(no_interp, idw, nn, las_merged)

    data_splitter(create_split, visualise_split, DATA_SPLIT_FOLDER_TO_VISUALIZE)
    
    get_png(FOLDER_LIST, NODATA, dtype)

    if len(HEIGHT_FOLDER_PATH) != 0:
        visualize_empty_height_data(HEIGHT_FOLDER_PATH, GEOJSON_EMPTY_TILES_PATH, NODATA_PNG)

    # end time
    end = time.time()
    print(f"Finished all tasks in {int((end - start)/60)} min, {(end - start)%60} sec.")