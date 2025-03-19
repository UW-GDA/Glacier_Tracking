import argparse
import time
import itertools
import numpy as np
import xarray as xr
import shapely
from tqdm import tqdm
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
import multiprocessing as mp
from glob import glob
from scipy.linalg import lstsq
from datetime import datetime as dt

from utils import get_auxiliary_angles


def get_glims_bounds(glimps_gdf, crs, bounds):
    # Clipping GLIMS polygons to Area of Interest
    glimps_gdf_aoi = glimps_gdf.clip(bounds)
    glimps_baker = glimps_gdf_aoi.to_crs(crs)
    glimps_baker['area'] = glimps_baker.area
    glimps_baker.sort_values('area', ascending=False, inplace=True)
    
    # Only consider the top 15 largest glaciers
    glimps_baker = glimps_baker.iloc[:15].buffer(30)
    return glimps_baker

def get_timesteps(offset_filenames):
    '''
    Get the timesteps from the offset filenames
    
    Parameters:
    ----------
    offset_filenames: list
        List of offset filenames
        
    Returns:
    -------
    timesteps: list
        Sorted list of timesteps
    '''
    timesteps = []
    for file in offset_filenames:
        timesteps.append(file[:-3].split('-')[0])
        timesteps.append(file[:-3].split('-')[1])

    timesteps = sorted(list(set(timesteps)))
    # print(len(timesteps), '22')
    # timesteps.pop(timesteps.index('20200326'))
    # print(len(timesteps), '33')
    return timesteps

def create_combined_dataset(off_filenames, offset_datapath, glimps_bounds, var):
    """
    Create a combined dataset of all the offset files

    Parameters:
        off_filenames (_type_): _description_
        offset_datapath (_type_): _description_
        glimps_bounds (_type_): _description_

    Returns:
        _type_: _description_
    """
    xr_files = []
    for i, file in tqdm(enumerate(off_filenames[:])):
        t1, t2 = file[:-3].split('/')[-1].split('-')
        fn = glob(f'{offset_datapath}/*/{t1}-{t2}*')[0]
        disp_da = xr.open_dataset(fn)
        disp_da.rio.write_crs(disp_da['spatial_ref'].attrs['spatial_ref'], inplace=True)
        disp_da = disp_da.rio.clip(glimps_bounds.geometry, disp_da.rio.crs)
        # xr_files[file] = disp_da
        xr_files.append(disp_da[var])
    
    xr_files = xr.concat(xr_files, pd.Index(off_filenames, name='filename'))
    return xr_files


def get_designmat(offset_filenames, timesteps, idx_i, comb_offsets, deltaT):
    A = np.zeros((len(offset_filenames), len(timesteps)-1))
    b = np.zeros((len(offset_filenames), 1))
    
    # strt = time.time()
    for i, file in enumerate(offset_filenames):
        t1, t2 = file[:-3].split('/')[-1].split('-')
        # fn = glob(f'{offset_datapath}/ascending/*/{t1}-{t2}*')[0]
        # disp_da = xr_files.sel() #xr_files[file]
        disp = comb_offsets.sel(filename=file, x=idx_i[1], y=idx_i[0]).values
        b[i] = disp
        # A[i, timesteps.index(t1)] = 1
        A[i, timesteps.index(t1):timesteps.index(t2)] = 1
    
    A = (A.T * deltaT).T
    A = A[~np.isnan(b[:,0]), :]
    b = b[~np.isnan(b[:,0])]
    
    return A, b

def compute_single_sbas(offset_filenames, timesteps, idx_i_list, comb_offsets, deltaT, i):
    '''
    Compute the SBAS inversion for a chunk of pixels
    '''
    # print(len(timesteps), i)
    V_nev = np.zeros((len(idx_i_list), len(timesteps)-1))
    for ix, idx_i in enumerate(idx_i_list):
        A, b = get_designmat(offset_filenames, timesteps, idx_i, comb_offsets, deltaT)
        
        v_skip_idx = []
        if (A.sum(0)==0).any():
            v_skip_idx = np.where(A.sum(0)==0)[0][0]
            A = np.delete(A, (v_skip_idx), axis=1)
            # continue
        
        # Checking if overdetermined or underdetermined
        if A.shape[0] > A.shape[1]:
            # v_skip_idxs.append(v_skip_idx)
            
            assert A.shape[0] == b.shape[0], 'Shape mismatch'
            
            v_nev, res, rank, s = lstsq(A, b, lapack_driver='gelsy', check_finite=False)
            v_nev = np.insert(v_nev, (v_skip_idx), np.nan)
        else:
            print('Underdetermined')
            v_nev = np.zeros((len(timesteps)-1))
            v_nev[:] = np.nan
        try:
            V_nev[ix] = v_nev
        except:
            # print(v_nev.shape, V_nev.shape, len(timesteps), v_skip_idx, ix, idx_i)
            V_nev[ix] = v_nev
    
    print(i)
    
    return V_nev

def compute_sbas(offset_filenames, test_aoi_dt, comb_offsets, var, offset_path, chunk_size=1000):
    '''
    Compute the SBAS inversion for all the pixels using multiprocessing
    
    Parameters:
    ----------
    offset_filenames: list
        List of offset filenames
    test_aoi_dt: xarray.Dataset
        Dataset containing the area of interest
    comb_offsets: xarray.Dataset
        Dataset containing the combined offsets for all the pixels
    var: str
        Variable to invert
    offset_path: str
        Path to the offset files
    chunk_size: int
        Number of pixels to process in a single iteration
    
    Returns:
    -------
    xr_files_out: xarray.Dataset
        Dataset containing the SBAS inversion results
    '''
    
    deltaT = [(dt.strptime(fn[:-3].split('-')[1], '%Y%m%d') - dt.strptime(fn[:-3].split('-')[0], '%Y%m%d')).days for fn in offset_filenames]
    timesteps = get_timesteps(offset_filenames)
    
    # creating empty DataArray to store the results
    xr_files_out = xr.zeros_like(comb_offsets)
    xr_files_out[:] = np.nan
    timesteps_out = [f'{timesteps[i]}-{timesteps[i+1]}.nc' for i in range(len(timesteps)-1)]
    xr_files_out = xr_files_out.reindex(filename=timesteps_out)

    idx_valid = np.argwhere((~np.isnan(test_aoi_dt[var])).values)
    idx_valid[:,0] = test_aoi_dt.y.values[idx_valid[:,0]]
    idx_valid[:,1] = test_aoi_dt.x.values[idx_valid[:,1]]

    num_cores = min(mp.cpu_count(), 18)
    chunk_inputs = [(offset_filenames, timesteps, idx_valid[j:j+chunk_size], comb_offsets, deltaT, j) for j in range(0, len(idx_valid)-chunk_size, chunk_size)]
    if len(idx_valid)%chunk_size != 0:
        chunk_inputs.append((offset_filenames, timesteps, idx_valid[-(len(idx_valid)%chunk_size):], comb_offsets, deltaT, len(idx_valid)-(len(idx_valid)%chunk_size)))
    
    with mp.Pool(num_cores) as pool:
        results = pool.starmap(compute_single_sbas, chunk_inputs)
    
    results = list(itertools.chain.from_iterable(results))
    
    strt = time.time()
    for i, idx_i in enumerate(idx_valid[:]):
        xr_files_out.loc[{'x': idx_i[1], 'y': idx_i[0]}] = results[i]
    print(time.time()-strt)
    
    xr_files_out.attrs['timesteps_out'] = timesteps_out
    
    return xr_files_out

def main():
    
    parser = argparse.ArgumentParser(description='Inversion of 3D displacements using SBAS')
    # arg.add_argument('--glims', type=str, help='Path to GLIMS shapefile')
    parser.add_argument('--offsets', type=str, help='Path to offset files')
    parser.add_argument('--var', type=str, help='Variable to invert [Options: "az" and "slrng"]')
    # parser.add_argument('--output', type=str, help='Path to output directory')
    args = parser.parse_args()    
    
    var = f'D_{args.var}'
    glimps_gdf = gpd.read_file('/mnt/Backups/ayushg12/GDA_proj/shapefiles/glims_polygons.shp')
    lats = np.array([-121.96, 48.713,-121.6869, 48.713,-121.6869, 48.8539,-121.96, 48.8539,-121.96, 48.713])[::2]
    lons = np.array([-121.96, 48.713,-121.6869, 48.713,-121.6869, 48.8539,-121.96, 48.8539,-121.96, 48.713])[1::2]
    bounds = np.array([lats.min(), lons.min(), lats.max(), lons.max()])
    
    offset_files = glob(f'{args.offsets}/deltaT_*/2021*.nc')
    offset_files = [fl for fl in offset_files if fl.split('/')[-2] not in ['deltaT_6', 'deltaT_18', 'deltaT_30']]
    offset_filenames = sorted([file.split('/')[-1] for file in offset_files])
    
    test_dt = xr.open_dataset(offset_files[0])
    test_dt.rio.write_crs(test_dt['spatial_ref'].attrs['spatial_ref'], inplace=True)

    glimps_baker = get_glims_bounds(glimps_gdf, test_dt.rio.crs, bounds)

    # Getting indexes for all valid pixels
    test_aoi_dt = test_dt.rio.clip(glimps_baker.geometry, test_dt.rio.crs)
    
    comb_offsets = create_combined_dataset(offset_filenames, args.offsets, glimps_baker, var)
    
    offset_sbas_out = compute_sbas(offset_filenames, test_aoi_dt, comb_offsets, var, args.offsets)
    
    offset_sbas_out.to_netcdf(f'{args.offsets}/sbas_inversion_{var}.nc')


if __name__=='__main__':
    main()
