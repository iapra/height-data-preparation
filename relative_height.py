import math
import numpy as np
import csv
import laspy
import os
import shapely.wkt
from shapely.geometry import Point, Polygon
import math
import numpy as np
import pandas as pd
import geopandas as gpd

def get_point_gdf(las_path):
    las_files = [os.path.join(las_path,l) for l in os.listdir(las_path) if l.endswith('.las')]
    lidar_points = []
    print('Reading .las files...')
    for las_file in las_files:
        with laspy.open(las_file) as f:
            las = f.read()
            building_points = las.points[las.classification == 6] # filter for building_pionts
            lidar_points.extend([Point(xyz) for xyz in zip(building_points.x,building_points.y,building_points.z)]) # append points
    print('Converting to GeoDataFrame...')
    df=pd.DataFrame(lidar_points,columns=['geometry'])
    gdf_lidar_points=gpd.GeoDataFrame(df,crs='epsg:25832',geometry='geometry')
    gdf_lidar_points['point']=gdf_lidar_points['geometry'] # duplicate geometry column to retain it in sjoin operation (sjoin drops active geometry column in right table)
    return gdf_lidar_points

def get_poly_gdf(POLY_CSV):
    df = pd.read_csv(POLY_CSV,sep=',')
    df['geom_3d'] = df['geom_3d'].apply(shapely.wkt.loads)
    gdf_polygons = gpd.GeoDataFrame(df,crs='epsg:25832',geometry='geom_3d')
    return gdf_polygons

def plane_equation_fit (points):
    '''
    Computes the plane equation of the plane best fitting given points.
    Output:
    [a,b,c,d] = coefficients of the plane equation, ax + by + cz = d 
    '''
    # G is the points barycenter (centered coordinates)
    G = points.sum(axis=0) / points.shape[0]
    # run SVD
    u, s, vh = np.linalg.svd(points - G)
    # unitary normal vector: gives a, b, c for ax + by + cz + d = 0
    u_norm = vh[2, :]
    a,b,c = u_norm
    # Get D (a * Xcentered_point + b * Ycentered_point + c * Zcentered_point which equals d)
    d = np.dot(u_norm, G)
    coef_equation = [a,b,c,d]
    return coef_equation

def isAbove_plane(point, equation_coef):
    '''
    Outputs if a point is above (True) or under a plane (False)
    Input: 
    point = point coordinates in array form [x,y,z]
    equation_coef = coefficients of the plane equation, in the form [a,b,c,d]
    '''
    a, b, c, d = equation_coef[0], equation_coef[1], equation_coef[2], -equation_coef[3]
    x,y,z = point[0], point[1], point[2]
    # IF THE ACTUAL Z VALUE IS HIGHER THAN THE EXPECTED ONE, 
    # THE POINT IS ABOVE THE ROOF-PLANE
    expected_z = -(a*x + b*y +d) /c
    if z >= expected_z:
        return True
    else: 
        return False

def shortest_distance(p, equation_coef):
    '''
    Computes shortest/minimum distance from point p to plane
    Input: 
    p = point coordinates in array form [x,y,z]
    equation_coef = coefficients of the plane equation, in the form [a,b,c,d]
    '''
    a, b, c, d = equation_coef[0], equation_coef[1], equation_coef[2], -equation_coef[3]
    x, y, z = p[0], p[1], p[2]
    numerator = abs((a * x + b * y + c * z + d))
    denum = (math.sqrt(a * a + b * b + c * c))
    #print(numerator/denum)
    isAbove_result = isAbove_plane(p, equation_coef)
    if isAbove_result == True:
        return numerator/denum
    if isAbove_result == False:
        return -(numerator/denum)

def vertical_distance(point, equation_coef):
    '''
    Computes vertical distance between a point and a plane given its equation
    Input:
    point = point coordinates in array form [x,y,z]
    equation_coef = coefficients of the plane equation, in the form [a,b,c,d]
    '''
    z = point[2]
    a, b, c, d = equation_coef[0], equation_coef[1], equation_coef[2], -equation_coef[3]
    x,y,z = point[0], point[1], point[2]
    # VERTICAL DISTANCE IS THE DIFFERENCE BEWTEEN Z AND EXPECTED-Z
    expected_z = -(a*x + b*y +d) /c
    vertical_distance = z - expected_z
    return vertical_distance

def write_new_las(input_las_folder, dict_points, classification):
    '''
    Writes one or numerous new LAS files based on LAS file inputs, 
    changing only their Z values according to the dictionary given as input.
    Input: 
    folder_las = directory containing one or numerous raw las files
    dict_points =   keys are the point id (corresponding to the order in which the points were read from the LAS file), 
                    values are the new Z values
    classification = LAS classification of interest (e.g. enter 6 for building points)
    '''
    files = os.listdir(input_las_folder)
    id = 0
    for f in files:
        if f[-4:] == ".las":
            las_fp = input_las_folder + "/" + f
            if os.path.exists(f'{las_fp[:-4]}_new.las'):
                print(f'New las for {las_fp} already exists. It will be replaced.')
                os.remove(f'{las_fp[:-4]}_new.las')
            else: 
                print("Creating new LAS file for ", f)
            point_array_x, point_array_y, point_array_z, classification_arr = [], [], [], []
            with laspy.open(las_fp) as fh:
                # classification = 6
                # print(f"Classification {classification}")
                las = fh.read()
                scales, offsets = las.header.scales, las.header.offsets
                # Create a new header
                header_new = laspy.LasHeader(point_format=1, version="1.2")
                header_new.offsets = las.header.offsets
                header_new.scales = las.header.scales               
                # header = las.header
                for p in range(len(las.points)):
                    if (las.points[p].raw_classification == classification):   # Ground == 2; Buildings == 6
                        x = (float(las.points[p].X * scales[0]) + offsets[0])
                        point_array_x.append(x)
                        y = (float(las.points[p].Y * scales[1]) + offsets[1])
                        point_array_y.append(y)
                        try:
                            z = dict_points[id]
                        except:
                            z = -9999
                        point_array_z.append(z)
                        classification_arr.append(classification)
                        id += 1
            pointrecord = laspy.ScaleAwarePointRecord.zeros(len(point_array_x), header=header_new)
            pointrecord.x = point_array_x
            pointrecord.y = point_array_y
            pointrecord.z = point_array_z
            pointrecord.classification = classification_arr
            
            # pointrecord.classification = classification
            with laspy.open(f'{las_fp[:-4]}_new.las', mode="w", header=header_new) as writer:
                writer.write_points(pointrecord)
            print(f"New file written in {las_fp[:-4]}_new.las")

def get_new_las_files(distance_type, input_las_folder, gdf_polys_to_points, classification):
    '''
    Obtains new las files giving relative height between LiDAR points and LOD model roof polygons
    Input:
    distance_type = 'vertical' or 'min', str
    input_las_folder = path to folder containing las file(s) input, str
    gdf_polys_to_points = geodataframe containing points and their underlying roof polygons,
    gdf
    classification = las classification of interest, int
    '''
    dict_poly = {}
    dict_points = {}
    # NEW Z COMPUTATION
    print('Computing relative distances...')
    for _, row in gdf_polys_to_points.iterrows():
        dist = -9999
        # RETRIEVE INFO IN GDF
        point_3d = shapely.wkt.loads(str(row['point']))
        poly_id = int(row['poly_id'])
        point_id = int(row['index_point'])
        # DISTANCE TO ROOF PLANE
        p = list(point_3d.coords)
        point = [p[0][0], p[0][1], p[0][2]]
        # WE STORE PLANE-EQUATIONS IN A DICT SO IT IS COMPUTED ONLY ONCE PER POLYGON
        if poly_id in dict_poly:
            coef_equation = dict_poly[poly_id]
        # WE COMPUTE PLANE EQUATION IF IT HAD NOT BEEN DONE YET
        else:
            poly_3d = shapely.wkt.loads(str(row['geom_3d']))
            points = np.array(poly_3d.exterior.coords)
            coef_equation = plane_equation_fit(points)
            dict_poly[poly_id] = coef_equation
        # MINIMUM OR VERTICAL DISTANCE
        if distance_type == 'min':
            dist = shortest_distance(point, coef_equation)
        if distance_type == 'vertical':
            dist = vertical_distance(point, coef_equation)

        dict_points[point_id] = dist
    
    # WRITE NEW .LAS WITH NEW Z VALUES
    write_new_las(input_las_folder, dict_points, classification)

    return
