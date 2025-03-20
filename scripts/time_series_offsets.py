# Import required packages
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from tqdm import tqdm
import xarray as xr

import logging
root_logger = logging.getLogger()
# root_logger.setLevel('WARNING')
root_logger.setLevel(logging.WARNING)
import geopandas as gpd
from compute_offsets import compute_offsets
import glob


def select_pairs(scene_list, min_temp_bline, max_temp_bline):
    '''
    Function to select pairs of scenes based on temporal baseline
    
    Parameters
    ----------
    scene_list: list
        List of scenes to select pairs from
    min_temp_bline: int
        Minimum temporal baseline in days
    max_temp_bline: int
        Maximum temporal baseline in days
    
    Returns
    -------
    date_pairs: list
        List of date pairs
    scene_pairs: list
        List of scene pairs
    '''
    scene_dates = {}
    for scene in scene_list:
        date = scene.split('_')[-5][:8]
        scene_dates[date] = scene

    date_pairs, scene_pairs = [], []
    for date1 in scene_dates:
        for date2 in scene_dates:
            if datetime.strptime(date2, '%Y%m%d')-datetime.strptime(date1, '%Y%m%d') < timedelta(days=max_temp_bline) and not datetime.strptime(date2, '%Y%m%d')-datetime.strptime(date1, '%Y%m%d') < timedelta(days=min_temp_bline) and not date1 >= date2 :
                date_pairs.append([date1, date2])
                scene_pairs.append([scene_dates[date1], scene_dates[date2]])
    
    print(f'number of pairs: {len(date_pairs)}')
    
    return date_pairs, scene_pairs


if __name__=="__main__":
    
    parser = argparse.ArgumentParser(description='Inversion of 3D displacements using SBAS')
    parser.add_argument('--glims_path', type=str, help='Path to GLIMS shapefile')
    parser.add_argument('--aoi_wkt', type=str, help='WKT string of Area of Interest')
    parser.add_argument('--data_path', type=str, help='Path to OPERA data')
    args = parser.parse_args()

    proj_path = '../'
    os.makedirs(f'{proj_path}/offsets', exist_ok=True)
    
    # Getting scene lists
    scene_list_asc = sorted(glob.glob(f'{args.data_path}/*T064-135608-IW1*'))
    scene_list_des = sorted(glob.glob(f'{args.data_path}/*T013-026559-IW1*'))
    
    # Reading GLIMPS shapefiles
    glimps_gdf = gpd.read_file(args.glims_path)
    lats = np.array([float(pnt.split(' ')[0]) for pnt in args.aoi_wkt[9:-2].split(',')])
    lons = np.array([float(pnt.split(' ')[1]) for pnt in args.aoi_wkt[9:-2].split(',')])
    bounds = np.array([lats.min(), lons.min(), lats.max(), lons.max()])
    glimps_gdf_aoi = glimps_gdf.clip(bounds)
    
    # Computing offsets for ascending pairs
    print("Computing offsets for ascending pairs")
    date_pairs_asc, scene_pairs_asc = select_pairs(scene_list_asc, 0, 37)
    os.makedirs(f'{proj_path}/offsets/ascending', exist_ok=True)
    for dates, pairs in tqdm(zip(date_pairs_asc, scene_pairs_asc)):
        deltaT = (datetime.strptime(dates[1], '%Y%m%d') - datetime.strptime(dates[0], '%Y%m%d')).days
        if os.path.exists(f'{proj_path}/offsets/ascending/deltaT_{deltaT}/{dates[0]}-{dates[1]}.nc'):
            continue
        
        offset_mltlk, _ = compute_offsets(pairs, glimps_gdf_aoi, 3)
        os.makedirs(f'{proj_path}/offsets/ascending/deltaT_{deltaT}', exist_ok=True)
        offset_mltlk.to_netcdf(f'{proj_path}/offsets/ascending/deltaT_{deltaT}/{dates[0]}-{dates[1]}.nc')
    
    # Computing offsets for descending pairs
    print("Computing offsets for descending pairs")
    date_pairs_des, scene_pairs_des = select_pairs(scene_list_des, 0, 37)
    os.makedirs(f'{proj_path}/offsets/descending', exist_ok=True)
    for dates, pairs in tqdm(zip(date_pairs_des, scene_pairs_des)):
        deltaT = (datetime.strptime(dates[1], '%Y%m%d') - datetime.strptime(dates[0], '%Y%m%d')).days
        if os.path.exists(f'{proj_path}/offsets/descending/deltaT_{deltaT}/{dates[0]}-{dates[1]}.nc'):
            continue
        
        offset_mltlk, _ = compute_offsets(pairs, glimps_gdf_aoi, 3)
        os.makedirs(f'{proj_path}/offsets/descending/deltaT_{deltaT}', exist_ok=True)
        offset_mltlk.to_netcdf(f'{proj_path}/offsets/descending/deltaT_{deltaT}/{dates[0]}-{dates[1]}.nc')
        
    print("Done!")
