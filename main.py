import time
import numpy as np

from height_grids_main import height_grids_main
from relative_height_main import relative_height_main
from data_splitter import data_splitter
from utils import get_png
from config import LAS_INPUT_FOLDER, DATA_SPLIT_FOLDER, \
    DIR_ABS_HEIGHT_GEOTIFF, DIR_REL_HEIGHT_GEOTIFF

# 1 # 
### UPDATE CONFIG FILE

######################## RELATIVE HEIGHT ################################
# 2 #
### INDICATE THE LAS CLASSIFICATION YOU WANT
classification = 6
### INDICATE WHAT YOU NEED WITH TRUE/FALSE
get_sql_cmd = False
get_intermediate_new_z = False
get_new_las_from_txt = False
get_new_las_directly = False

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
FOLDER_LIST = []
NODATA = -9999
dtype = np.uint8

########################################################################

if __name__ == "__main__":
    # starting time
    start = time.time()

    relative_height_main(get_sql_cmd, get_intermediate_new_z, \
        get_new_las_from_txt, classification, get_new_las_directly)

    height_grids_main(no_interp, idw, nn, las_merged)

    data_splitter(create_split, visualise_split, DATA_SPLIT_FOLDER_TO_VISUALIZE)
    
    get_png(FOLDER_LIST, NODATA, dtype)

    # end time
    end = time.time()
    print(f"Finished all tasks in {int((end - start)/60)} min, {(end - start)%60} sec.")