# Import required packages
import os
import numpy as np
from datetime import datetime, timedelta
import xarray as xr
import logging
root_logger = logging.getLogger()
# root_logger.setLevel('WARNING')
root_logger.setLevel(logging.WARNING)
import isce
from components.contrib.geo_autoRIFT.autoRIFT import autoRIFT_ISCE
from scipy.interpolate import interpn
os.environ['AUTORIFT'] = '/mnt/Backups/gbrench/sw/insar_tools/isce2/src/isce2/contrib/geo_autoRIFT'
os.environ['PYTHONPATH'] = os.environ['AUTORIFT']
os.environ['PATH'] += f":{os.environ['AUTORIFT']}"


def run_autoRIFT(pair, aoi_gdf, skip_x=3, min_azm_chip=16, max_azm_chip=64,
                 preproc_filter_width=21, mpflag=10, search_limit_x=4):

    print(f'opening pair: {pair[0]}-{pair[1]}')
    
    ds1 = xr.open_dataset(pair[0], group='data', engine='h5netcdf')
    crs = ds1.projection.attrs['spatial_ref']
    aoi_gdf_proj = aoi_gdf.to_crs(crs)
    aoi_bounds = aoi_gdf_proj.total_bounds
    aoi_bounds[:2] -= 8000
    aoi_bounds[2:] += 8000
    
    x_spacing = np.abs(ds1.x_coordinates[1:].values - ds1.x_coordinates[:-1].values).mean()
    y_spacing = np.abs(ds1.y_coordinates[1:].values - ds1.y_coordinates[:-1].values).mean()
    ratio_yx = np.ceil(y_spacing/x_spacing)
    ds1 = ds1.rio.write_crs(crs)
    ds1 = ds1.rio.clip_box(*aoi_bounds)
    slc1 = ds1.VV.values

    ds2 = xr.open_dataset(pair[1], group='data', engine='h5netcdf')
    ds2 = ds2.rio.write_crs(crs)
    ds2 = ds2.rio.clip_box(*aoi_bounds)
    slc2 = ds2.VV.values

    I1 = np.abs(slc1)
    I1[np.isnan(I1)] = 0
    I2 = np.abs(slc2)
    I2[np.isnan(I2)] = 0
    
    obj = autoRIFT_ISCE()
    obj.configure()
    obj.MultiThread = mpflag

    # rotate to take advantage of chip size scaling
    I1 = np.rot90(I1)
    I2 = np.rot90(I2)
    # scale range chip size to get nearly square chip in cartesian coordinates
    obj.ScaleChipSizeY = ratio_yx

    obj.I1 = I1
    obj.I2 = I2

    obj.SkipSampleX = skip_x
    obj.SkipSampleY = int(np.ceil(skip_x*ratio_yx))

    # Kernel sizes to use for correlation
    obj.ChipSizeMinX = min_azm_chip
    obj.ChipSizeMaxX = max_azm_chip
    obj.ChipSize0X = min_azm_chip
    # oversample ratio, balancing precision and performance for different chip sizes
    obj.OverSampleRatio = {obj.ChipSize0X:32,obj.ChipSize0X*2:64,obj.ChipSize0X*4:128}

    # generate grid
    m,n = obj.I1.shape
    xGrid = np.arange(obj.SkipSampleX+10,n-obj.SkipSampleX,obj.SkipSampleX)
    yGrid = np.arange(obj.SkipSampleY+10,m-obj.SkipSampleY,obj.SkipSampleY)
    nd = xGrid.__len__()
    md = yGrid.__len__()
    obj.xGrid = np.int32(np.dot(np.ones((md,1)),np.reshape(xGrid,(1,xGrid.__len__()))))
    obj.yGrid = np.int32(np.dot(np.reshape(yGrid,(yGrid.__len__(),1)),np.ones((1,nd))))
    noDataMask = np.invert(np.logical_and(obj.I1[:, xGrid-1][yGrid-1, ] > 0, obj.I2[:, xGrid-1][yGrid-1, ] > 0))
    
    # set y search limit based on ratio_yx
    search_limit_y = int(np.ceil(search_limit_y*ratio_yx))
    # set search limits
    obj.SearchLimitX = np.full_like(obj.xGrid, search_limit_x)
    obj.SearchLimitY = np.full_like(obj.xGrid, search_limit_y)

    # set search limit and offsets in nodata areas
    obj.SearchLimitX = obj.SearchLimitX * np.logical_not(noDataMask)
    obj.SearchLimitY = obj.SearchLimitY * np.logical_not(noDataMask)
    obj.Dx0 = obj.Dx0 * np.logical_not(noDataMask)
    obj.Dy0 = obj.Dy0 * np.logical_not(noDataMask)
    obj.Dx0[noDataMask] = 0
    obj.Dy0[noDataMask] = 0
    obj.NoDataMask = noDataMask

    print("preprocessing")
    obj.WallisFilterWidth = preproc_filter_width
    obj.preprocess_filt_hps()
    obj.uniform_data_type()

    print("starting autoRIFT")
    obj.runAutorift()
    print("autoRIFT complete")

    # return outputs to original orientation
    tmpDx = np.rot90(obj.Dx, axes=(1, 0))
    obj.Dx = np.rot90(obj.Dy, axes=(1, 0))
    obj.Dy = tmpDx
    tmpxGrid = np.rot90(obj.xGrid, axes=(1, 0))
    obj.xGrid = np.rot90(obj.yGrid, axes=(1, 0))
    obj.yGrid = tmpxGrid
    obj.InterpMask = np.rot90(obj.InterpMask, axes=(1, 0))
    obj.ChipSizeX = np.rot90(obj.ChipSizeX, axes=(1, 0))
    obj.NoDataMask = np.rot90(obj.NoDataMask, axes=(1, 0))
    obj.bounds = aoi_bounds
    obj.x_coords_geo = ds1.x_coordinates.values
    obj.y_coords_geo = ds1.y_coordinates.values
    obj.crs = crs

    # convert displacement to m
    # obj.Dx_m = obj.Dx * 2.3
    # obj.Dy_m = obj.Dy * 12.1
    obj.Dx_m = obj.Dx * x_spacing
    obj.Dy_m = obj.Dy * y_spacing
        
    return obj


def compute_offsets(scene_pairs, aoi_gdf, skip_x=4, search_limit_x=4):
    '''
    Compute offsets between two scenes using autoRIFT
    
    Parameters
    ----------
    scene_pairs : list
        List of two scene pairs to compare
    aoi_gdf : geopandas.GeoDataFrame
        Area of interest polygon
    skip_x : int
        Skip sampling distance in x direction for autoRIFT
    
    Returns
    -------
    offsets_multilooked : xarray.Dataset
        Multilooked offsets
    offsets : xarray.Dataset
        Full resolution offsets
    '''
    obj = run_autoRIFT([scene_pairs[0], scene_pairs[1]], aoi_gdf, skip_x=skip_x, search_limit_x=search_limit_x)

    x_coords = np.flip(obj.xGrid[0, :])
    y_coords = obj.yGrid[:, 0]

    # Create a mesh grid for the slc dimensions
    x_coords_new, y_coords_new = np.meshgrid(
        np.arange(obj.I2.shape[0]),
        np.arange(obj.I2.shape[1])
    )

    print("Starting interpolation")
    # Perform bilinear interpolation using scipy.interpolate.interpn
    Dx_full = interpn((y_coords, x_coords), obj.Dx_m, (y_coords_new, x_coords_new), method="linear", bounds_error=False)
    Dy_full = interpn((y_coords, x_coords), obj.Dy_m, (y_coords_new, x_coords_new), method="linear", bounds_error=False)
    chipsize_full = interpn((y_coords, x_coords), obj.ChipSizeX, (y_coords_new, x_coords_new), method="linear", bounds_error=False)
    interpmask_full = interpn((y_coords, x_coords), obj.InterpMask, (y_coords_new, x_coords_new), method="linear", bounds_error=False)
    print("Interpolation complete")

    x_coords_proj = obj.x_coords_geo
    y_coords_proj = obj.y_coords_geo

    offsets = xr.Dataset(
        {"D_slrng": (("y", "x"), Dx_full), "D_az": (("y", "x"), Dy_full),
         "chip_size_azm":(['y', 'x'], chipsize_full), "interpMask":(['y', 'x'], interpmask_full)},
        {"x": x_coords_proj, "y": y_coords_proj},
        attrs={"description": "Offset tracking results using autoRIFT on OPERA products",
            "pair_1": scene_pairs[0].split('/')[-1],
            "pair_2": scene_pairs[1].split('/')[-1],
            "Displacement_units": "meters"},
    )

    offsets = offsets.rio.write_crs(obj.crs)

    # # multilook to match igram
    offsets_multilooked = offsets.coarsen(x=2, boundary='trim').mean()
    
    return offsets_multilooked, offsets

