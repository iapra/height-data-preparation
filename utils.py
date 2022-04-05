from osgeo import gdal
from gdalconst import GA_ReadOnly
import laspy
import numpy as np

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

        # [ulx_deg, uly_deg, np.float128(lrx_deg), np.float128(lry_deg)]
        # [np.float128(ulx_deg), np.float128(uly_deg), np.float128(lrx_deg), np.float128(lry_deg)]

def get_array(las_fp):
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