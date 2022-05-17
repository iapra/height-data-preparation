import os
import geopandas as gpd
import pickle as pkl

from config import LAS_INPUT_FOLDER, POLY_CSV

from relative_height import get_point_gdf, get_poly_gdf, get_new_las_files

def relative_height_main(distance_type, classification, get_new_las):
    '''
    Central function to obtain LAS file(s) with relative height
    between LiDAR  point cloud and LOD roof polygons
    Input:
    distance_type = 'vertical' or 'min', str
    classification = las classification of interest, int
    get_new_las = bool
    '''
    if get_new_las == True:
        # WE STORE THE GDF IN PKL FOR FASTER RE-LOAD IF NECESSARY
        pkl_file = f"{LAS_INPUT_FOLDER}/gdf_polys_to_points.pkl"
        if os.path.exists(pkl_file):
            with open(pkl_file, 'rb') as f:
                [gdf_polys_to_points] = pkl.load(f)
        elif not os.path.exists(pkl_file):
            # WE STORE POINTS AND ROOF POLYGONS IN GEODATAFRAMEs
            gdf_lidar_points = get_point_gdf(LAS_INPUT_FOLDER)
            gdf_polygons = get_poly_gdf(POLY_CSV)
            # WE COUPLE EACH POINT TO A POLYGON
            gdf_polys_to_points = gpd.tools.sjoin(gdf_polygons, gdf_lidar_points, how='inner',rsuffix='point',lsuffix='polygon')
            with open(pkl_file, 'wb') as f:
                pkl.dump([gdf_polys_to_points], f)
        get_new_las_files(distance_type, LAS_INPUT_FOLDER, gdf_polys_to_points, classification)
