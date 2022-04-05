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
    Function to find distance from point p to plane
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
    Function to find distance from point p to plane
    Input: 
    p = point coordinates in array form [x,y,z]
    equation_coef = coefficients of the plane equation, in the form [a,b,c,d]
    '''
    a, b, c, d = equation_coef[0], equation_coef[1], equation_coef[2], -equation_coef[3]
    x, y, z = p[0], p[1], p[2]
    numerator = abs((a * x + b * y + c * z + d))
    denum = (math.sqrt(a * a + b * b + c * c))
    #print(numerator/denum)
    return numerator/denum

def txt_SQL_cmd(folder_las, file_out, table, classification):
    '''
    This function create an SQL command to import LAS points to a DB
    input: folder_las = directory containing one or numerous raw las files
    file_out = directory for the txt file containing the SQL command
    table = table name in which to enter the points
    classification = LAS classification of interest (e.g. enter 6 for building points)
    '''
    files = os.listdir(folder_las)
    array_txt = []
    id = 0
    for f in files:
        if f[-4:] == ".las":
            las_fp = folder_las + "/" + f
            print(f"Adding points from {f}, starting with id {id}")
            with laspy.open(las_fp) as fh:
                las = fh.read()
                scales, offsets = las.header.scales, las.header.offsets
                for p in range(len(las.points)):
                    if (las.points[p].classification == classification):   # Ground == 2; Buildings == 6
                        x = (float(las.points[p].X * scales[0]) + offsets[0])
                        y = (float(las.points[p].Y * scales[1]) + offsets[1])
                        z = (float(las.points[p].Z * scales[2]) + offsets[2])
                        value_strg = f"({id}, ST_GeomFromEWKT('SRID=25832;POINT({x} {y} {z})'))"
                        array_txt.append(value_strg)
                        id += 1
    print(len(array_txt))
    with open(file_out, 'w') as f:
        f.write(f"INSERT INTO {table} (point_id, geom) \nVALUES ")
        for values in array_txt[:-1]:
            f.write(f"{values},\n")
        f.write(f"{array_txt[-1]}")
        f.close()

def write_new_las(input_las_folder, dict_points, classification):
    '''
    This function write one or numerous new LAS files based on LAS file inputs, 
    changing only their Z values according to the dictionary given as input.
    input: folder_las = directory containing one or numerous raw las files
    dict_points =   keys are the point id (corresponding to the order in which the points were read from the LAS file), 
                    values are the new Z values
    classification = LAS classification of interest (e.g. enter 6 for building points)

    NB --- folder_las and classification have to be similar to the ones used for the txt_SQL_cmd
    so that the point ids are correct
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
                        
def get_new_z(file_out_z, POLY_CSV, POINTS_CSV, txt_output_bool):
    dict_poly = {}
    dict_points = {}
    # READ THE ROOF-POLYGONS CSV FILE, CONTAINING POLY_ID, WKT GEOM_2D, WKT GEOM_3D
    with open(POLY_CSV) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                # print(f'Column names are {", ".join(row)}')
                line_count += 1
                continue
            else:
                poly_id, geom_3d, geom_2d = row
                poly_2d = shapely.wkt.loads(geom_2d)
                poly_3d = shapely.wkt.loads(geom_3d)
                # WE COMPUTE THE PLANE EQUATION FITTING BEST THE POLYGON POINTS
                points = np.array(poly_3d.exterior.coords)
                coef_equation = plane_equation_fit(points)
                # WE STORE EVERYTHING IN A DICT
                dict_poly[poly_id] = (poly_2d, poly_3d, coef_equation)
                line_count += 1
        print(f'Processed {line_count} polygons.')
        
    with open(POINTS_CSV) as csv_file2:
        csv_reader = csv.reader(csv_file2, delimiter=',')
        line_count = 0
        for point in csv_reader:
            # if line_count == 31046:
            #     return
            if line_count == 0:
                # print(f'Column names are {", ".join(point)}')
                line_count += 1
                continue
            else:
                if line_count%100000 == 0:
                    print(f"{line_count} points done...")
                point_id, point_geom_3d, point_geom_2d = point
                point_2d = shapely.wkt.loads(point_geom_2d)
                point_3d = shapely.wkt.loads(point_geom_3d)
                # RETRIEVE UNDERLYING POLYGON
                min_dist = -9999
                for poly_id in dict_poly:
                    poly_2d, poly_3d, coef_equation = dict_poly[poly_id]
                    if poly_2d.contains(point_2d):
                        # COMPUTE DISTANCE TO PLANE
                        min_dist = point_3d.distance(poly_3d)
                        p = list(point_3d.coords)
                        point = [p[0][0], p[0][1], p[0][2]]
                        min_dist = shortest_distance(point, coef_equation)
                        # THE DISTANCE IS NEGATIVE IF THE POINT IS UNDER THE PLANE
                        isAbove_result = isAbove_plane(point, coef_equation)
                        if isAbove_result == True:
                            if txt_output_bool == True:
                                with open(file_out_z, 'a') as f:
                                    f.write(f"{point_id}, {min_dist}\n")
                            else:
                                dict_points[point_id] = min_dist
                        if isAbove_result == False:
                            if txt_output_bool == True:
                                with open(file_out_z, 'a') as f:
                                    f.write(f"{point_id}, {-min_dist}\n")
                            else:
                                dict_points[point_id] = -min_dist
                        break
                    else:
                        continue
                # THE DISTANCE IS SET TO -9999 (NODATA) IF THERE IS NO PLANE UNDER THE POINT
                if min_dist == -9999:
                    if txt_output_bool == True:
                        with open(file_out_z, 'a') as f:
                            f.write(f"{point_id}, {min_dist}\n")
                    else:
                        dict_points[point_id] = min_dist
                line_count += 1

        print(f'Processed {line_count} points.')
        # if txt_output_bool == False:
        #     assert line_count == len(dict_points)
        return dict_points

def get_new_las(input_las_folder, file_out_z, classification):
    # READ NEW Z VALUES FROM TXT AND CREATE DICT
    dict_points = {}
    with open(file_out_z) as f:
        for line in f:
            content = line.strip().split(', ')
            id, z = int(content[0]), float(content[1])
            dict_points[id] = z
    print("Dict read ! Writing las file(s)...")

    # WRITE NEW .LAS WITH NEW Z VALUES
    write_new_las(input_las_folder, dict_points, classification)