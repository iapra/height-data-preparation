import os
from osgeo import gdal
import numpy as np
import scipy.spatial
import laspy
import os
import rasterio
import scipy.ndimage
import scipy.interpolate
import imageio
from PIL import Image
import rasterio as rio

from config import IDW_power, IDW_radius, IMG_SIZE

from utils import get_array

def get_merge_pipeline(las_folder_path, las_merged, in_srs, out_srs):
    string_input = ""
    for img_name in os.listdir(las_folder_path):
        if img_name[-4:] != ".las":
            continue
        else:
            img_fp = las_folder_path + "/" + img_name
            # if len(string_input) == 0:
            string_input = string_input + '\
                {"type":"readers.las", \
                "filename":"' + img_fp + '"},\n'
    command = '[\
            ' + string_input + '\n \
    { "type": "filters.merge" },\n \
    { "type":"filters.reprojection",\n\
        "in_srs":"EPSG:' + str(in_srs) + '",\n\
        "out_srs":"EPSG:' + str(out_srs) + '"\n\
        },\
    { \
        "type":"writers.las", \
        "filename":"' + las_merged + '" \
    } \
        ]'
    print(command)
    return command

def get_crop_pipeline(xmin, xmax, ymin, ymax, las_merged, las_folder_path, geotiff_in):
    command = '[\
    "' + las_merged + '",\
    { "type":"filters.crop",\
        "bounds":"([' + str(xmin) + ',' + str(xmax) + '], [' + str(ymin) + ',' + str(ymax) + '])",\
        "a_srs": "EPSG:4326"\
    },\
        { \
        "type":"writers.las", \
        "filename":"' + las_folder_path + '/buildings_cropped/' + str(geotiff_in) + '.las" \
    } \
     ]'
    return command

def get_no_interp_pipeline(las_cropped, gtif_out_no_interp, cellSize, origin, size):
    command = '[\
    "' +  las_cropped + '",\
    { "type":"writers.gdal",\
        "filename":"' + gtif_out_no_interp + '",\
            "output_type":"max",\
                "gdaldriver":"GTiff",\
                    "radius":0.14,\
                        "power":1,\
                            "resolution":' + str(cellSize) + ',\
                                "origin_x":' + str(origin[0]) + ',\
                                    "origin_y":' + str(origin[1]) + ',\
                                        "height":' + str(size) + ',\
                                            "width":' + str(size) + ',\
                            "nodata": -9999\
    }\
    ]'
    return command

def get_idw_pipeline(las_cropped, gtif_out_idw, IDW_power, IDW_radius, cellSize, origin, size):
    command = '[\
    "' +  las_cropped + '",\n\
    { "type":"writers.gdal",\n\
        "filename":"' + gtif_out_idw + '",\n\
            "output_type":"idw",\n\
                "gdaldriver":"GTiff",\n\
                    "radius": %s,\n\
                        "power": %s,\n\
                            "resolution": %s,\n\
                                "origin_x": %s,\n\
                                    "origin_y": %s,\n\
                                        "height": %s,\n\
                                            "width": %s,\n\
                            "nodata": -9999,\n\
                                "where": "!(Z == -9999)"\
    }\
    ]' % (IDW_radius, IDW_power, cellSize, origin[0], origin[1], size, size)
    return command


def nn_interpolation(img_fp, bbox_img, p, gtif_out_nn):
    """
    Function that writes the output raster with nearest neighbour interpolation
     
    Input:
        p: the array of the input points (in 3D)
    Output:
        returns the value of the area
    """  
    cellSize = (bbox_img[1][0]- bbox_img[0][0])/IMG_SIZE
    #-- to speed up the nearest neighbour search, we use a kd-tree
    points = np.array(p)

    # First, we define the coordinates of the centers of the cells of the output grid (X,Y)
    ll, ur = bbox_img
    arrayX = np.linspace(ll[0]+ cellSize/2 , ur[0]+ cellSize/2, num=IMG_SIZE)
    arrayY = np.linspace(ll[1]+ cellSize/2, ur[1]+ cellSize/2, num=IMG_SIZE)
    x, y = np.meshgrid(arrayX, arrayY)
    centers = np.array(list(zip(x.flatten(),y.flatten())))

    # Now, we compute the distance from cell-center to the sample point
    X_tree = scipy.spatial.KDTree(points[:,:2])
    distance, indices = X_tree.query(centers)

    # If the cell is outside the convex-hull it will have "no-data" value,
    # otherwise it will get the z value of the closest sample point.
    hull = scipy.spatial.Delaunay(points[:,:2])
    z = []
    i = 0
    for index in indices:
        if hull.find_simplex(centers[i])>=0 and distance[i]<=1.5:    # the point is in the convex-hull
            z.append(points[index][2])
        else :                                                       # the point is outside the convex-hull
            z.append(-9999) 
        i += 1
    z = np.reshape(np.array(z), newshape = (len(arrayY), len(arrayX)))
    # z2 = z[::-1,:]

    # TIF FILE OUT
    with rio.open(img_fp) as src:
        ras_meta = src.profile
        ras_meta.update(count=1,
        dtype="float64",
        nodata=-9999)
    with rio.open(gtif_out_nn, 'w', **ras_meta) as dst:
        dst.write(z, 1)
    print("File written to", gtif_out_nn)
    
def get_nn_interpolation(img_fp, bbox_img, las_cropped, gtif_out_nn):
    if os.path.exists(gtif_out_nn):
        print("Already done!")
        return
    else:
        array_bdg = get_array(las_cropped)
        if len(array_bdg) == 0:
            get_empty_gtif(img_fp, gtif_out_nn)
        else:
            nn_interpolation(img_fp, bbox_img, array_bdg, gtif_out_nn)


def get_empty_gtif (fp_img, gtif_out):
    with rio.open(fp_img) as src:
        ras_meta = src.profile
    with rio.open(gtif_out, 'w', **ras_meta) as dst:
        empty_array = np.zeros(shape=(IMG_SIZE,IMG_SIZE))
        empty_array[empty_array==0] = np.nan
        dst.write(empty_array, 1)
    print("Empty file written to", gtif_out)

def interpolation_grids(img_fp, las_cropped, interpolation_folder, \
                        img_name, bbox_img, origin, \
                        no_interp, idw, nn):
    cellSize = (bbox_img[1][0]- bbox_img[0][0])/IMG_SIZE
    # GRIDS PATH OUT
    gtif_out_no = f"{interpolation_folder}/no_interpolation/{img_name[:-4]}.tif"
    gtif_out_idw = f"{interpolation_folder}/idw_interpolation/{img_name[:-4]}.tif"
    # asc_out_nn = f"{interpolation_folder}/nn_interpolation/{img_name[:-4]}.asc"
    gtif_out_nn = f"{interpolation_folder}/nn_interpolation/{img_name[:-4]}.tif"
    if no_interp == True:
        if not os.path.exists(f"{interpolation_folder}/no_interpolation"):
            os.mkdir(f"{interpolation_folder}/no_interpolation")
        if not os.path.exists(f"{interpolation_folder}/no_interpolation/temp"):
            os.mkdir(f"{interpolation_folder}/no_interpolation/temp")
        if os.path.exists(gtif_out_no):
            print("No-interpolation already exists!")
        # # CREATE AND RUN PDAL NN-PIPELINE
        if not os.path.exists(gtif_out_no):
            json_out_no_interp = f"{interpolation_folder}/no_interpolation/temp/{img_name[:-4]}.json"
            cmd1 = get_no_interp_pipeline(las_cropped, gtif_out_no, cellSize, origin, size=IMG_SIZE)
            with open (json_out_no_interp, "w") as outfile_no_interp:
                outfile_no_interp.write(cmd1)
            try:
                os.system(f"pdal pipeline {json_out_no_interp}")
            except:
                get_empty_gtif(img_fp, gtif_out_no)
            
    
    if idw == True:
        if not os.path.exists(f"{interpolation_folder}/idw_interpolation"):
            os.mkdir(f"{interpolation_folder}/idw_interpolation")
        if not os.path.exists(f"{interpolation_folder}/idw_interpolation/temp"):
            os.mkdir(f"{interpolation_folder}/idw_interpolation/temp")
        if os.path.exists(gtif_out_idw):
            print("IDW-interpolation already exists !")
        if not os.path.exists(gtif_out_idw):
            json_out_idw = f"{interpolation_folder}/idw_interpolation/temp/{img_name[:-4]}.json"
            # CREATE AND RUN PDAL IDW-PIPELINE
            cmd2 = get_idw_pipeline(las_cropped, gtif_out_idw, \
                                    IDW_power, IDW_radius, \
                                    cellSize, origin, size=IMG_SIZE)
            with open (json_out_idw, "w") as outfile_idw:
                outfile_idw.write(cmd2)
            try:
                os.system(f"pdal pipeline {json_out_idw}")
            except:
                get_empty_gtif(img_fp, gtif_out_idw)
    
    if nn == True:
        if not os.path.exists(f"{interpolation_folder}/nn_interpolation"):
            os.mkdir(f"{interpolation_folder}/nn_interpolation")
        # NN interpolation
        get_nn_interpolation(img_fp, bbox_img, las_cropped, gtif_out_nn)
    
    return



        