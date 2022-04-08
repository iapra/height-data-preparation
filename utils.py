from osgeo import gdal
from gdalconst import GA_ReadOnly
import laspy
import numpy as np
import csv
import os
from shapely.ops import cascaded_union
from shapely.geometry import Polygon
import geopandas as gpd

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
