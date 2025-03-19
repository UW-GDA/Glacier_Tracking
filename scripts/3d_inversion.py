import argparse
import time
import itertools
import numpy as np
import xarray as xr
import shapely
from tqdm import tqdm
import multiprocessing as mp
from glob import glob
from scipy.linalg import lstsq
from datetime import datetime as dt

from utils import get_matA_3D
from compute_sbas import get_timesteps
from utils import get_auxiliary_angles


def get_timesteps_clean(asc_Daz, des_Daz):

    deltaT_gap = (des_Daz.time.values[0] - asc_Daz.time.values[0])
    asc_timesteps = [val + deltaT_gap for val in asc_Daz.time.values]

    timesteps_all = sorted(list(set(list(asc_timesteps) + list(des_Daz.time.values))))
    timesteps_intersect = sorted([val for val in timesteps_all if val in asc_timesteps and val in list(des_Daz.time.values)])
    timesteps_exclude = [val for val in timesteps_all if val not in timesteps_intersect]
    
    excl_list_asc, excl_list_des = [], []
    for excl in timesteps_exclude:
        if excl in asc_timesteps:
            ascT = excl - deltaT_gap
            # ascT = dt.strftime(ascT, '%Y%m%d')
            excl_list_asc.append(ascT)
        elif excl in des_Daz.time.values:
            # desT = dt.strftime(excl, '%Y%m%d')
            excl_list_des.append(excl)
            
    return excl_list_asc, excl_list_des, timesteps_intersect

def clean_data(Daz, Dslrng, excl_list):
    # Cleaning ascending data
    for excl in excl_list:
        # Getting index of the time step
        idx = np.where(Daz.time.values==excl)[0][0]
        # Merging idx index with previous index
        Daz.loc[{'time': Daz.time.values[idx-1]}] += Daz.loc[{'time': excl}]
        Dslrng.loc[{'time': Dslrng.time.values[idx-1]}] += Dslrng.loc[{'time': excl}]
        # Removing the time step
        Daz = Daz.drop_sel(time=Daz.time.values[idx-1])
        Dslrng = Dslrng.drop_sel(time=Dslrng.time.values[idx-1])
    
    return Daz, Dslrng

def get_timesteps_out(offset_path):
    asc_filenames = glob(f'{offset_path}/*/2021*.nc')
    offset_filenames = sorted([file.split('/')[-1] for file in asc_filenames])
    timesteps = [dt.strptime(val, '%Y%m%d') for val in get_timesteps(offset_filenames)]
    
    return timesteps

def read_offset_data(offset_path):
    timesteps = get_timesteps_out(offset_path)
    deltaT = np.array([(timesteps[i+1] - timesteps[i]).days  for i in range(len(timesteps)-1)])/365 # For converting to m/yr

    Daz = xr.open_rasterio(f'{offset_path}/sbas_inversion_D_az.nc')
    Daz = Daz.assign_coords(band=timesteps[:-1])
    Daz = Daz.rename({'band':'time'})
    Dslrng = xr.open_rasterio(f'{offset_path}/sbas_inversion_D_slrng.nc')
    Dslrng = Dslrng.assign_coords(band=timesteps[:-1])
    Dslrng = Dslrng.rename({'band':'time'})

    # Coverting descending displacements to velocity [m/yr]
    Daz = Daz/deltaT.reshape(-1, 1, 1)
    Dslrng = Dslrng/deltaT.reshape(-1, 1, 1)
    
    return Daz, Dslrng

def compute_3DVel_single_hermet(V_aoro, A):
    A1, A2 = A[::2], A[1::2]
    b1, b2 = V_aoro[::2], V_aoro[1::2]
    P = np.eye(4)
    # P = initP(np.diag(P), V_aoro)
    P1 = np.diag(np.diag(P)[::2]) # LOS
    P2 = np.diag(np.diag(P)[1::2]) # Azimuth
    for i in range(3):
        N = (A1.T@P1@A1) + (A2.T@P2@A2)
        U = (A1.T@P1@b1) + (A2.T@P2@b2)
        v_nev, res, rank, s = lstsq(N, U, lapack_driver='gelsy', check_finite=False)
        resi1 = (b1 - np.dot(A1, v_nev)) #.reshape(4)    
        sigma01 = (resi1.T@P1@resi1)/2
        
        resi2 = (b2 - np.dot(A2, v_nev)) #.reshape(4)
        sigma02 = (resi2.T@P2@resi2)/2
        if sigma02==0:
            break
        
        P2 = (sigma01/sigma02)*P2
        
        if np.abs(sigma01-sigma02)<1:
            break
    
    res = (V_aoro - np.dot(A, v_nev))
    return v_nev.reshape(3), res 

def compute_3DVel_single(V_aoro, A):
    
    v_nev, res, rank, s = lstsq(A, V_aoro, lapack_driver='gelsy', check_finite=False)
    return v_nev.reshape(3), res 

def unpacking_chunk(jj, idx_valid_list, A, vel_comb_asc, vel_comb_des):
    v_nev = np.empty((len(idx_valid_list), len(vel_comb_asc.time.values), 3), dtype=np.float32)
    v_nev.fill(np.nan)
    
    for ti in range(len(vel_comb_asc.time.values)):
        # vel_comb_ti = vel_comb.isel(time=ti)
        for ix, idx_i in enumerate(idx_valid_list):
            v_aoro_asc = vel_comb_asc.isel(time=ti).sel(x=idx_i[1], y=idx_i[0]).values.reshape(-1, 1)
            v_aoro_des = vel_comb_des.isel(time=ti).sel(x=idx_i[1], y=idx_i[0]).values.reshape(-1, 1)
            v_aoro = np.vstack([v_aoro_asc, v_aoro_des])
            # v_nev[ix, ti, :], _ = compute_3DVel_single_hermet(v_aoro, A)
            v_nev[ix, ti, :], _ = compute_3DVel_single(v_aoro, A)

    print(f'Chunk {jj} done')
    return v_nev

def compute_3DVel_mp(asc_Daz, asc_Dslrng, des_Daz, des_Dslrng, aux_data, chunk_size=1000):
    A = get_matA_3D(*aux_data)
    # timesteps = asc_Daz.time.values.shape[0]
    vel_comb_asc = xr.concat([asc_Dslrng, asc_Daz], dim='band')
    vel_comb_des = xr.concat([des_Dslrng, des_Daz], dim='band')

    mask = ~(np.isnan(asc_Daz.values).any(axis=0)|np.isnan(asc_Dslrng.values).any(axis=0)|np.isnan(des_Daz.values).any(axis=0)|np.isnan(des_Dslrng.values).any(axis=0))
    idx_valid = np.argwhere(mask)
    idx_valid[:,0] = asc_Daz.y.values[idx_valid[:,0]]
    idx_valid[:,1] = asc_Daz.x.values[idx_valid[:,1]]

    num_cores = min(mp.cpu_count(), 18)
    chunk_inputs = [(j, idx_valid[j:j+chunk_size], A, vel_comb_asc, vel_comb_des) for j in range(0, len(idx_valid)-chunk_size, chunk_size)]
    if len(idx_valid)%chunk_size != 0:
        chunk_inputs.append((len(idx_valid)-(len(idx_valid)%chunk_size), idx_valid[-(len(idx_valid)%chunk_size):], A, vel_comb_asc, vel_comb_des))

    with mp.Pool(num_cores) as pool:
        results = pool.starmap(unpacking_chunk, chunk_inputs)

    results = list(itertools.chain.from_iterable(results))
    print("Inversion completed")
    
    # creating empty DataArray to store the results
    vel3d_out_SN = xr.zeros_like(asc_Daz)
    vel3d_out_SN[:] = np.nan
    vel3d_comb = xr.merge([vel3d_out_SN.rename('velocity_SN'), vel3d_out_SN.copy().rename('velocity_WE'), vel3d_out_SN.copy().rename('velocity_V')])
    
    for i, idx_i in tqdm(enumerate(idx_valid[:])):
        vel3d_comb['velocity_SN'].loc[{'x': idx_i[1], 'y': idx_i[0]}] = results[i][:, 0]
        vel3d_comb['velocity_WE'].loc[{'x': idx_i[1], 'y': idx_i[0]}] = results[i][:, 1]
        vel3d_comb['velocity_V'].loc[{'x': idx_i[1], 'y': idx_i[0]}] = results[i][:, 2]

    return vel3d_comb



def main():
    
    parser = argparse.ArgumentParser(description='Inversion of 3D displacements using SBAS')
    parser.add_argument('--datapath', type=str, help='Path to data directory')
    args = parser.parse_args()
    
    # Shapely polygon of Area of Interest
    baker_bounds = shapely.from_wkt('POLYGON((-121.96 48.713,-121.6869 48.713,-121.6869 48.8539,-121.96 48.8539,-121.96 48.713))')

    # Get incidence angles for ascending and descending passes
    asc_incid, asc_azim  = get_auxiliary_angles(baker_bounds, 64)
    des_incid, des_azim = get_auxiliary_angles(baker_bounds, 13)
    asc_azim = 360 + asc_azim
    des_azim = 360 + des_azim
    aux_data = [np.deg2rad(asc_azim), np.deg2rad(asc_incid), np.deg2rad(des_azim), np.deg2rad(des_incid)]
    
    asc_offset_path = f'{args.datapath}/ascending'
    des_offset_path = f'{args.datapath}/descending'

    asc_Daz, asc_Dslrng = read_offset_data(asc_offset_path)
    des_Daz, des_Dslrng = read_offset_data(des_offset_path)
    
    # Aligning the time steps of ascending and descending data
    exclude_list_asc, exclude_list_des, timesteps_interesect = get_timesteps_clean(asc_Daz, des_Daz)
    asc_Daz, asc_Dslrng = clean_data(asc_Daz, asc_Dslrng, exclude_list_asc)
    des_Daz, des_Dslrng = clean_data(des_Daz, des_Dslrng, exclude_list_des)
    
    print("Computing 3D velocities")
    vel3d_comb = compute_3DVel_mp(asc_Daz, asc_Dslrng, des_Daz, des_Dslrng, aux_data)
    vel3d_comb.to_netcdf(f'{args.datapath}/3d_vel.nc')


if __name__ == "__main__":
    main()
