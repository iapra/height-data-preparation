# Relative & Absolute Height data grids

## Description
A repository to prepare height data for training a DL model. Absolute height data is derived from LiDAR point clouds directly. Relative height data is derived from subtracting LiDAR points to LOD model available. 

## Container
The pip package requirements are listed in the "requirements.txt".  
Requirements can be setup through the container built with the provided Dockerfile. The latter includes python libraries, GDAL-3.1.4 and PDAL-2.3.0.  

Setup, once docker is installed on your machine:  
*$ cd path/to/project-folder*  
*$ docker build -t myname/myproject:tag .*  
The whole directory can be mounted using -v flag:  
*$ docker run -it -v $(pwd):/app/height-data \           
myname/myproject:tag*  

Once the container is run, scripts can be run from the shell:    
*$ python height-data/main.py*

## Data
Data is not provided in this repository. Required data include digital orthophotos (DOPs), LiDAR point cloud (LAS format) and LOD model in case of relative height computation. The DOPs can be acquired through WMS request, using the corresponding script. The paths to input data should be updated in the config file.

## Content

* main.py  
Gathers all the main functions with True/False parameters to update.
* DOP_grid.py  
Generates a DOP dataset from WMS request.
The photos are building centered. Buildings' centres are extracted from Open Street Map.
* relative_height_main.py  
Computes the difference between a LiDAR point cloud and its corresponding LOD model.
New LAS files are generated with relative Z values.  
--> Uses functions defined in *relative_height.py*.
* height_grids_main.py  
Create interpolation grids of the input LAS files, based on the boundary of the DOPs images.  
--> Uses functions defined in *height_grids.py*.
* data_splitter.py  
Outputs folder with 3 txt files: train.txt, val.txt, test.txt containing image ids so that datasets do not overlap.
* config.py  
Indicates WMS address, paths and interpolation configurations.
* utils.py  
* config.py  
--> TO ADD
* data/  
--> TO ADD

### Config file

Add a config file with following content:  

#### # BASE DIRECTORY IN CONTAINER  
DIR_BASE = *str*  
##### # WMS passwords and address  
username = *str*  
password = *str*  
wms_adress = *str*  
#### # Interpolations parameters  
IDW_radius = *int*  
IDW_power = *int*  
IMG_SIZE = *int*      # square image  
#### # SQL COMMANDS TO IMPORT LAS POINT DATA TO DATABASE  
FILE_OUT_SQL = *str*  
TABLE = *str*  
#### # CSV OUTPUT/INPUT for normalised/relative height data  
POINTS_CSV = *str*  
POLY_CSV = *str*  
#### # RELATIVE Z OUTPUT TXT-FILE  
FILE_OUT_Z = *str*  
#### # PATHS FOR .TIF, .LAS AND INTERPOLATION  
INTERPOLATION_FOLDER = *str*  
DIR_REL_HEIGHT_GEOTIFF = INTERPOLATION_FOLDER + *str* 
DIR_ABS_HEIGHT_GEOTIFF = INTERPOLATION_FOLDER + *str* 
LAS_INPUT_FOLDER = *str*  
DIR_IMAGES_GEOTIFF = *str*  
#### # TRAIN, VAL, TEST SPLITS  
DATA_SPLIT_FOLDER = *str*       # main name from which will be made several split_xxx folders  
NUMBER_OF_SPLITS = *int*  
PERCENTAGE_TRAIN, PERCENTAGE_VAL, PERCENTAGE_TEST = *float, float, float - values in [0-1]*  