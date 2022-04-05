import os
from tkinter import image_names
from owslib.wms import WebMapService
from owslib.util import Authentication
import geopandas as gpd
from shapely.geometry import Point, Polygon
#create config.py and add it to .gitignore
from config_superstructures import username, password, PATH_OUT, DIR_IMAGES_GEOTIFF 
import numpy as np
import gdal
from gdalconst import GA_ReadOnly

def get_bbox(input_path):
    data = gdal.Open(input_path, GA_ReadOnly)
    geoTransform = data.GetGeoTransform()
    minx = geoTransform[0]
    maxy = geoTransform[3]
    maxx = minx + geoTransform[1] * data.RasterXSize
    miny = maxy + geoTransform[5] * data.RasterYSize
    # print([minx, miny, maxx, maxy])
    return data.RasterXSize, data.RasterYSize, [minx, miny, maxx, maxy]
    
def get_aerial_img_Geodaten_Bayern(username, password, pointlist_lon_lat, img_width, img_height, filepath):
    """
    Sets up a wms connection to geodaten bayern and saves an aerial image from 
    a list of shapely points
    
    Parameters
    ----------
    username : string
    
    password : string
    
    pointlist_lon_lat : list of shapely.geometry.Point with [LON, LAT]
    
    img_size : total number of pixels (max. 4000x4000)
    
    filepath : string
        Filepath for image storage.
    
    """
    # set up authentication
    user, pw = username, password
    authenticate = Authentication()
    
    #connect to service
    wms = WebMapService('https://geoservices.bayern.de/wms/v2/ogc_dop20.cgi?',
                        version='1.1.1',
                        xml=None, 
                        username=user,
                        password=pw, 
                        parse_remote_metadata=False, 
                        timeout=30,
                        headers=None, 
                        auth=authenticate)
    
    # define bounding box using the min and max values of x and y    
    gs = gpd.GeoSeries(pointlist_lon_lat)
    # bottom left
    blc = list([gs.bounds.minx.min(),gs.bounds.miny.min()]) 
    # top right
    trc = list([gs.bounds.maxx.max(),gs.bounds.maxy.max()]) 
    
    # get the tile and save as image    
    tile = wms.getmap(layers=['by_dop20c'],
                styles=['default'],
                srs='EPSG:4326',
                bbox = (blc[0],blc[1],trc[0],trc[1]),
                size=(img_width, img_height),
                format='image/png',
                transparent=True)   
    
    out = open(filepath, 'wb')
    out.write(tile.read())
    out.close()

    return

for img_name in os.listdir(DIR_IMAGES_GEOTIFF):
    filepath = PATH_OUT + f"/images_DOP_roof_centered_png/{img_name}"
    if img_name[-4:] != ".tif":
        # print(f"I continue ! {img_name}")
        continue
    if os.path.exists(filepath):
        # print(f"already exists ! {img_name}")
        continue
    else:
        print(f"doing ! {img_name}")
        fp = DIR_IMAGES_GEOTIFF + "/" + img_name
        img_width, img_height, [xmin, ymin, xmax, ymax] = get_bbox(fp)
        res_diff = 10/20 # google res/DOP resolution
        # res_diff=1
        
        p1 = Point([xmin, ymin])
        p2 = Point([xmin, ymax])
        p3 = Point([xmax, ymax])
        p4 = Point([xmax, ymin])

        pointlist_lon_lat = [p1, p2, p3, p4]

        get_aerial_img_Geodaten_Bayern(username, password,
        pointlist_lon_lat, int(img_width*res_diff), int(img_height*res_diff), filepath)