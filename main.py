import time

from height_grids_main import height_grids_main
from relative_height_main import relative_height_main
from data_splitter import data_splitter
from config import LAS_INPUT_FOLDER, DATA_SPLIT_FOLDER

# starting time
start = time.time()

# 1 # 
### UPDATE CONFIG FILE

######################## RELATIVE HEIGHT ################################
# 2 #
### INDICATE THE LAS CLASSIFICATION YOU WANT
classification = 6      # Refers to LAS classification rules

# 3 #
### INDICATE WHAT YOU NEED WITH TRUE/FALSE
relative_height_main(get_sql_cmd = False, \
    get_intermediate_new_z = False, get_new_las_from_txt = False,\
    classification = classification, get_new_las_directly = False)

######################### INTERPOLATION #################################
# 2 #
### INDICATE A SINGLE MERGED POINT CLOUD OR A PATH TO CREATE IT
las_merged = f"{LAS_INPUT_FOLDER}/merged.las"

# 3 #
### ENTER WHAT INTERPOLATION YOU WANT WITH TRUE/FALSE
height_grids_main(no_interp = False, idw = False, nn =True, \
    las_merged = las_merged)

######################### DATA-SPLIT ###################################
# 2 #
### Enter the directory of split folder to visualise if necessary
DATA_SPLIT_FOLDER_TO_VISUALIZE = DATA_SPLIT_FOLDER + "_3"

# 3 #
# INDICATE WHAT YOU NEED WITH TRUE/FALSE
data_splitter(create_split = False, visualise_split = False, \
    DATA_SPLIT_FOLDER_TO_VISUALIZE = DATA_SPLIT_FOLDER_TO_VISUALIZE)
########################################################################

# end time
end = time.time()
print(f"Runtime of the program was {int((end - start)/60)} min, {(end - start)%60} sec.")