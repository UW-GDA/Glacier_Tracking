import ee
import math
import datetime as dt
import geopandas as gpd
import numpy as np
import json

# Authenticate and initialize GEE
try:
    ee.Initialize(project='ee-ayushg12')
except Exception as e:
    ee.Authenticate()
    ee.Initialize(project='ee-ayushg12')
    
    
def get_auxiliary_angles(bounds, relative_orbit):
    """
    Function to compute the mean incidence angle over a given region of interest (ROI) using Sentinel-1 SAR data
    
    Parameters
    ----------
    bounds: shapely.geometry.Polygon
        Shapely polygon representing the region of interest
    relative_orbit: int
        Relative orbit number of the Sentinel-1 data
        
    Returns
    -------
    angle_value: float
        Mean incidence angle over the ROI
    """
    # Convert Shapely polygon to GeoJSON
    geojson_poly = json.loads(gpd.GeoSeries([bounds]).to_json())["features"][0]["geometry"]

    roi_np = np.array(geojson_poly['coordinates'])[:,:-1]
    
    # Convert GeoJSON to an Earth Engine Geometry
    roi = ee.Geometry.Polygon(roi_np.tolist(), proj='EPSG:4326')

    # Load Sentinel-1 SAR collection and filter by date, region, and polarization
    sentinel1 = (ee.ImageCollection("COPERNICUS/S1_GRD")
                .filterBounds(roi)
                .filterDate("2020-01-01", "2021-01-01")  # Adjust date range
                .filter(ee.Filter.eq("relativeOrbitNumber_stop", relative_orbit))
                .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))  # Ensure VV is available
                .select("angle"))  # Select the incidence angle band

    # Compute mean incidence angle over ROI
    incidence_angle = sentinel1.median().reduceRegion(
        reducer=ee.Reducer.median(),
        geometry=roi,
        scale=10,  # Adjust resolution if needed
        maxPixels=1e13
    )

    # Get and print the result
    angle_value = incidence_angle.getInfo().get("angle")
    
    heading_angle = np.median(sentinel1.aggregate_array("platformHeading").getInfo())

    
    return angle_value, heading_angle


def get_Va_mat(phi):
    return np.array([np.cos(phi), np.sin(phi), 0])


def get_Vr_mat(phi, incid):
    return np.array([(-np.sin(phi)*np.sin(incid)), (np.sin(incid)*np.cos(phi)), -np.cos(incid)])


def get_matA_3D(phi_asc, incid_asc, phi_des, incid_des):
    A = np.zeros((4, 3))
    A[0] = get_Vr_mat(phi_asc, incid_asc)
    A[1] = get_Va_mat(phi_asc)
    A[2] = get_Vr_mat(phi_des, incid_des)
    A[3] = get_Va_mat(phi_des)
    return A
    
def get_matB(velR_asc, velA_asc, velR_des, velA_des):
    b = np.zeros((4, 1))
    b[0] = velR_asc
    b[1] = velA_asc
    b[2] = velR_des
    b[3] = velA_des
    return b

# def get_files_dict(files):
#     file_dict = {}
#     for file in files:
#         d1, d2 = file[:-3].split('/')[-1].split('-')[0], file[:-3].split('/')[-1].split('-')[1]
#         if d1 not in file_dict.keys():
#             file_dict[d1] = {}
        
#         file_dict[d1][d2] = file

#     return file_dict
