from osgeo import gdal
from gdalconst import GA_ReadOnly
import laspy
import numpy as np
import csv
import os
import shutil
import tifffile as tiff
from shapely.ops import cascaded_union
from shapely.geometry import Polygon
import geopandas as gpd
import cv2
from PIL import Image

from config import DIR_IMAGES_GEOTIFF, IMG_SIZE

def bbox(points):
    """ Defines the two oposite corners of the bounding box 
    Input : points (array form)
    Returns two points (corners of the bbox)
    """
    lx = min(point[0] for point in points)
    ly = min(point[1] for point in points)
    ux = max(point[0] for point in points)
    uy = max(point[1] for point in points)
    p1 = (lx, ly)
    p2 = (ux, uy)
    bbox = [p1, p2]
    return bbox

def get_img_bbox(input_path):
    '''
    Obtains geographical information about .tif image

    Input:
    input_path = .tif image input path

    Output:
    RasterXSize = raster size X direction, 
    RasterYSize = raster size Y direction, 
    [minx, miny, maxx, maxy] = bouding-box projected coordinates,
    [ulx_deg, uly_deg, lrx_deg, lry_deg] = bouding-box coordinates in degrees
    '''
    ### GET PROJECTED COORDINATES
    data = gdal.Open(input_path, GA_ReadOnly)
    geoTransform = data.GetGeoTransform()
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1] * data.RasterXSize
    miny = maxy + geoTransform[5] * data.RasterYSize

    ### GET DEGREE COORDINATES
    wkt_srs = data.GetProjection()
    ulx, uly = gdal.ApplyGeoTransform(geoTransform, 0, 0)
    lrx, lry = gdal.ApplyGeoTransform(geoTransform, data.RasterXSize, data.RasterYSize)

    src_srs = gdal.osr.SpatialReference()
    src_srs.ImportFromWkt(wkt_srs)
    tar_srs = gdal.osr.SpatialReference()
    tar_srs.ImportFromEPSG(3857)

    # with recent versions of GDAL the axis order (x,y vs y,x) depends
    # on the projection. Force "x,y" with:
    src_srs.SetAxisMappingStrategy(gdal.osr.OAMS_TRADITIONAL_GIS_ORDER)
    tar_srs.SetAxisMappingStrategy(gdal.osr.OAMS_TRADITIONAL_GIS_ORDER)

    ct = gdal.osr.CoordinateTransformation(src_srs, tar_srs)

    ulx_deg, uly_deg, z = ct.TransformPoint(ulx, uly)
    lrx_deg, lry_deg, z = ct.TransformPoint(lrx, lry)
    # print(ulx_deg, uly_deg, lrx_deg, lry_deg )

    return data.RasterXSize, data.RasterYSize, [minx, miny, maxx, maxy], \
        [np.float128(ulx_deg), np.float128(uly_deg), np.float128(lrx_deg), np.float128(lry_deg)]

def get_array(las_fp):
    '''
    Function to obtain an array from .las input file
    Input:
    las_fp = path to las input file, str
    Output:
    array of 3D points [[x,y,z], [x,y,z], ...]
    '''
    array = []
    with laspy.open(las_fp) as fh:
        las = fh.read()
        scales, offsets = las.header.scales, las.header.offsets
        for p in range(len(las.points)):
            # if (las.points[p].classification == 6):
            x = (float(las.points[p].X * scales[0]) + offsets[0])
            y = (float(las.points[p].Y * scales[1]) + offsets[1])
            z = (float(las.points[p].Z * scales[2]) + offsets[2])
            arr = [x,y,z]
            array.append(arr)
    return array

def datasets_to_geojson(split_folder_path, images_dir):
    '''
    Outputs a geojson polygon per dataset (train, validation, test).
    Input:
    split_folder_path = path to split folder containing txt files with image_ids
    images_dir = geotif images directory, to associate image_ids to image boundaries
    '''
    train_txt, val_txt, test_txt = f"{split_folder_path}/train.txt", f"{split_folder_path}/val.txt", f"{split_folder_path}/test.txt"
    train_json, val_json, test_json = f"{split_folder_path}/train.json", f"{split_folder_path}/val.json", f"{split_folder_path}/test.json" 
    
    # We write train_images to geojson
    if not os.path.exists(train_json):
        _txt_to_geojson(train_txt, train_json, images_dir)

    # We write val_images to geojson
    if not os.path.exists(val_json):
        _txt_to_geojson(val_txt, val_json, images_dir)

    # We write test_images to geojson
    if not os.path.exists(test_json):
        _txt_to_geojson(test_txt, test_json, images_dir)

def _txt_to_geojson(txt_path, geojson_path, images_dir):
    '''
    Outputs a geojson per txt file contained in a split folder.
    Input:
    txt_path = path to txt file with image_ids
    geojson_path = output path for the geojson polygon
    images_dir = geotif images directory, to associate image_ids to iamge boundaries
    '''
    if not os.path.exists(geojson_path):
        polys = []
        with open(txt_path) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter='\n')
            # Each row of the txt is an image-id
            for image_id in csv_reader:
                image_path = os.path.join(images_dir, f"{str(image_id[0])}.tif")
                poly = _get_poly_from_img(image_path)
                polys.append(poly)
        train_merge = gpd.GeoSeries(cascaded_union(polys))
        train_merge.to_file(geojson_path, driver="GeoJSON") 

def _get_poly_from_img(img_path):
    '''
    Get image boundaries as polygon, given the .tif image path
    '''
    _, _, [minx, miny, maxx, maxy], _ = get_img_bbox(img_path)
    poly = Polygon([[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy]])
    return poly

def _las_point_number (las_path):
    '''
    Reads a .las file and returns number of points contained
    '''
    with laspy.open(las_path) as fh:
        las = fh.read()
        return (len(las.points))

def get_png(folder_list, NODATA, dtype):
    for FOLDER in folder_list:
        if os.path.exists(f"{FOLDER}_png"):
            shutil.rmtree(f"{FOLDER}_png")
            print("Removing existing png conversion.")

        if not os.path.exists(f"{FOLDER}_png"):
            os.mkdir (f"{FOLDER}_png")

        for infile in os.listdir(f"{FOLDER}"):
            if infile[-4:] == ".tif":
                fp = f"{FOLDER}/{infile}"
                outfile = infile.split('.')[0] + '.png'
                path_out = f"{FOLDER}_png/{outfile}"
                ## CONVERT FROM FLOAT64 TO UINT8
                if os.path.exists(path_out):
                    continue
                ## READ AND KEEP ONLY 3 BANDS FOR RGB
                try:
                    data = tiff.imread(fp)[:,:,:3]
                ## ELSE, ONLY ONE BAND
                except:
                    data = tiff.imread(fp)
                ## RETRIEVE MIN AND MAX VALUES FOR NORMALISING
                try:
                    min_value = min(data[data!=NODATA])
                except:
                    min_value = np.nan
                try:    
                    max_value = max(data[data!=NODATA]) - min_value
                except:
                    max_value = np.nan
                data = np.nan_to_num(data)
                if len(data[data==NODATA]) != 0:
                    try:
                        data[data==NODATA] = np.nan
                    except:
                        assert data.all() == 0
                        data_zeros = np.zeros((IMG_SIZE,IMG_SIZE))
                        tiff.imwrite(path_out, data_zeros.astype(dtype))
                        continue

                if min_value!=np.nan and max_value != 0:
                    data = ((data - min_value)/max_value)*255
                
                if min_value!=np.nan and max_value == 0:
                    data = ((data - min_value))*255

                tiff.imwrite(path_out, data.astype(dtype))

def visualize_empty_height_data(folder_path, geojson_path, no_data=0):
    """
    Outputs a geojson vector of all height images containing 
    no height data information, i.e. containing only no_data value

    Input: folder_path = path to the folder containing png images to evaluate, string
    geojson_path = output path for the geojson vector file
    no_data = value of no_data in the png images, int
    """
    empty_imgs = []
    height_imgs = os.listdir(folder_path)
    for img in height_imgs:
        if img[-4:] == ".png":
            height_path = os.path.join(folder_path, img)
            img_content = np.array(Image.open(height_path))
            if len(img_content[img_content != no_data]) == 0:
                image_path = f"{DIR_IMAGES_GEOTIFF}/{img[:-4]}.tif"
                poly = _get_poly_from_img(image_path)
                empty_imgs.append(poly)
            else:
                # The height image is not empty
                continue
        else:
            # This is not an image
            continue
    print(len(empty_imgs))
    empty_height_merged = gpd.GeoSeries(cascaded_union(empty_imgs))
    empty_height_merged.to_file(geojson_path, driver="GeoJSON") 
