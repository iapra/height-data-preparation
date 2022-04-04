## Relative & Absolute Height data grids

# Description
A repository to prepare height data for training a DL model. Absolute height data is derived from LiDAR point clouds directly. Relative height data is derived from subtracting LiDAR points to LOD model available. 

# Container
The pip package requirements are listed in the "requirements.txt".
Requirements can be setup through the container built with the provided Dockerfile. The latter includes python libraries, GDAL-3.1.4 and PDAL-2.3.0.

Setup, once docker is installed on your machine:
$ cd path/to/project-folder
$ docker build -t myname/myproject:tag .
The whole directory can be mounted using -v flag:
$ docker run -it -v $(pwd):/app/height-data \           
myname/myproject:tag

Once the container is run, scripts can be run from the shell:
$ python height-data/script.py

# Data
Data is not provided in this repository. Required data include digital orthophotos (DOPs), LiDAR point cloud (LAS format) and LOD model in case of relative height computation. The DOPs can be acquired through WMS request, using the corresponding script. The paths to input data should be updated in the config file.

# Content

.
 * DOP_grid
Generates a DOP dataset from WMS request.
The photos are building centered. Buildings' centres are extracted from Open Street Map.
 * relative_height_main
Computes the difference between a LiDAR point cloud and its corresponding LOD model.
New LAS file are generated with relative Z values.
 * interpolation_grid_main
Create interpolation grids of the input LAS files, based on the boundary of the DOPs.
 * main
 * config
Indicates WMS address, paths and interpolation configurations.
 * utils