import os
import numpy as np
import scipy.spatial
import os
import scipy.ndimage
import scipy.interpolate
import rasterio as rio
import math

from config import IDW_power, IDW_radius, IMG_SIZE

from utils import get_array, _las_point_number

def get_merge_pipeline(las_folder_path, las_merged, in_srs, out_srs):
    # First, we gather all input las files to merge
    string_input = ""
    for img_name in os.listdir(las_folder_path):
        if img_name[-4:] != ".las":
            continue
        else:
            img_fp = las_folder_path + "/" + img_name
            string_input = string_input + '\
                {"type":"readers.las", \
                "filename":"' + img_fp + '"},\n'
    # The string output is used as input for the pdal pipeline
    command = '[\
            ' + string_input + '\n \
    { "type": "filters.merge" },\n \
    { "type":"filters.reprojection",\n\
        "in_srs":"EPSG:' + str(in_srs) + '",\n\
        "out_srs":"EPSG:' + str(out_srs) + '"\n\
        },\n\
    {\n \
        "type":"writers.las",\n \
        "filename":"' + las_merged + '"\n \
    } \
        ]'
    return command

def get_crop_pipeline(xmin, xmax, ymin, ymax, las_merged, las_folder_path, geotiff_in):
    command = '[\n\
    "' + las_merged + '",\n\
    { "type":"filters.crop",\n\
        "bounds":"([' + str(xmin) + ',' + str(xmax) + '], [' + str(ymin) + ',' + str(ymax) + '])",\n\
        "a_srs": "EPSG:4326"\n\
    },\n\
        { \
        "type":"writers.las",\n \
        "filename":"' + las_folder_path + '/buildings_cropped/' + str(geotiff_in) + '.las"\n \
    }\n \
     ]'
    return command

def get_no_interp_pipeline(las_cropped, gtif_out_no_interp, cellSize, origin, size):
    # To get no interpolation, a radius of diagonal of a cellsize is used, and a power of 1
    radius = (cellSize/2)*math.sqrt(2)
    command = '[\
    "' +  las_cropped + '",\n\
    { "type":"writers.gdal",\n\
        "filename":"' + gtif_out_no_interp + '",\n\
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
    ]' % (radius, 1, cellSize, origin[0], origin[1], size, size)
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
    Function that writes the output raster with nearest neighbour interpolation of points
     
    Input:
    img_fp = file-path to the tif image, str
    bbox_img = ((lower-left-x, lower-left-y), (upper-right-x, upper-right-y))
    p = array of the input points in 3D, [[x,y,z],[x,y,z],...]
    gtif_out_nn = file-path for the output tif image

    Output:
    tif image for nearest-neighbor interpolation of the input 3D points
    """  
    cellSize = (bbox_img[1][0]- bbox_img[0][0])/IMG_SIZE
    points = np.array(p)

    # First, we define the coordinates of the centers of the cells of the output grid (X,Y)
    ll, ur = bbox_img
    arrayX = np.linspace(ll[0]+ cellSize/2 , ur[0]+ cellSize/2, num=IMG_SIZE)
    arrayY = np.linspace(ll[1]+ cellSize/2, ur[1]+ cellSize/2, num=IMG_SIZE)
    x, y = np.meshgrid(arrayX, arrayY)
    centers = np.array(list(zip(x.flatten(),y.flatten())))

    # To speed up the nearest neighbour search, we use a kd-tree
    # And compute the distances and indices of the point nearest to each cell-center
    X_tree = scipy.spatial.KDTree(points[:,:2])
    distance, indices = X_tree.query(centers)

    # If the cell is outside the convex-hull or too far from the cell-center, it will have "no-data" value,
    # otherwise it will get the z value of the closest sample point.
    hull = scipy.spatial.Delaunay(points[:,:2])
    z = []
    i = 0
    for index in indices:
        if hull.find_simplex(centers[i])>=0 and distance[i]<=1.5:    # the point is in the convex-hull & closer than 1.5 meter
            z.append(points[index][2])
        else :                                                       # the point is outside the convex-hull
            z.append(-9999) 
        i += 1
    assert len(z) == IMG_SIZE*IMG_SIZE
    z = np.reshape(np.array(z), newshape = (len(arrayY), len(arrayX)))
    z2 = z[::-1,:]

    # TIF FILE OUT
    with rio.open(img_fp) as src:
        ras_meta = src.profile
        ras_meta.update(count=1,
        dtype="float64",
        nodata=-9999)
    with rio.open(gtif_out_nn, 'w', **ras_meta) as dst:
        dst.write(z2, 1)
    print("File written to", gtif_out_nn)


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
    # GRIDS PATHS OUT
    gtif_out_no = f"{interpolation_folder}/no_interpolation/{img_name[:-4]}.tif"
    gtif_out_idw = f"{interpolation_folder}/idw_interpolation/{img_name[:-4]}.tif"
    gtif_out_nn = f"{interpolation_folder}/nn_interpolation/{img_name[:-4]}.tif"

    if no_interp == True:
        if not os.path.exists(f"{interpolation_folder}/no_interpolation"):
            os.mkdir(f"{interpolation_folder}/no_interpolation")
        if not os.path.exists(f"{interpolation_folder}/no_interpolation/temp"):
            os.mkdir(f"{interpolation_folder}/no_interpolation/temp")
        if os.path.exists(gtif_out_no):
            print("No-interpolation already exists!")
        # CREATE AND RUN PDAL NO-INTERPOLATION-PIPELINE
        if not os.path.exists(gtif_out_no):
            # IF THE POINT-CLOUD IS EMPTY, WRITE EMPTY TIF IMG
            if _las_point_number(las_cropped) == 0:
                get_empty_gtif(img_fp, gtif_out_no)
            # OTHERWISE, CREATE NO-INTERPOLATION
            else:
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
            # IF THE POINT-CLOUD IS EMPTY, WRITE EMPTY TIF IMG
            if _las_point_number(las_cropped) == 0:
                get_empty_gtif(img_fp, gtif_out_idw)
            # OTHERWISE, CREATE IDW-INTERPOLATION
            else:
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
        if os.path.exists(gtif_out_nn):
            print("NN-interpolation already exists !")
        if not os.path.exists(gtif_out_nn):
            # IF THE POINT-CLOUD IS EMPTY, WRITE EMPTY TIF IMG
            if _las_point_number(las_cropped) == 0:
                get_empty_gtif(img_fp, gtif_out_nn)
            # OTHERWISE, CREATE NN-INTERPOLATION
            else:
                try:
                    array_bdg = get_array(las_cropped)
                    nn_interpolation(img_fp, bbox_img, array_bdg, gtif_out_nn)
                except:
                    get_empty_gtif(img_fp, gtif_out_nn)
    
    return



        