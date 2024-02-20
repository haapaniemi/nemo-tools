#!/usr/bin/env python
"""
haapanie
02/02/2024

Extracting NEMO points from float trajectory
"""
import xarray as xr
import pandas as pd
import numpy as np
import argopy

if __name__ == '__main__':
    # Define path to data
    path = '/scratch/project_2003984/Veera/run_nemo'

    # Select time range (currently used for FLOAT data only)
    start_date = '2021-08-01'
    end_date = '2021-08-31'

    # Read in 6h grid T data
    exp = xr.open_dataset(path + '/CTRL_bfr/CTRL_exp_bfr_1d_20210401_20211001_grid_T_202108-202108.nc', engine="netcdf4") 

    # If already in a file:
    #     float = xr.open_dataset(path + '/results/float_6903710.nc')
    # Otherwise using argopy:
    float_no = 6903710
    dat = argopy.DataFetcher().float(float_no)
    float = dat.to_xarray()


    # Select time range
    fdata = float.to_dataframe()
    time_mask = (fdata['TIME'] >= start_date) & (fdata['TIME'] <= end_date)
    fdata = fdata[time_mask]

    # Drop duplicates --> only one time, lon, lat for each FLOAT CYCLE
    # (lat & lon same for each float cycle!)
    fdata = fdata.drop_duplicates('CYCLE_NUMBER')
    print('float data, dropped duplicates (dim = CYCLE_NUMBER):')
    print(fdata)
    fdata = fdata[['TIME', 'LONGITUDE', 'LATITUDE', 'TEMP', 'PRES']]
    print('float data reduced:')
    print(fdata)

    # Store float variables to data arrays
    flats = fdata['LATITUDE'].values
    flons = fdata['LONGITUDE'].values
    ftims = fdata['TIME'].values
    print('float time stamps:')
    print(ftims)

    # Create empty lists to store the x and y indices to
    xidx, yidx = [], []

    print('looping through float timestamps:')

    # Find the nearest indices for all of the time steps
    for i in range(len(ftims)):
        timestamp = ftims[i]
        print('........ current time stamp:')
        print(timestamp)

        # Select float lon & lat
        lon = flons[i]
        lat = flats[i]
        print('........ float lon, lat:')
        print(lon, lat)

        distance = np.sqrt((exp.nav_lon - lon)**2 + (exp.nav_lat - lat)**2)
        print('........ distance set shape (only used for findind nearest point):')
        print(distance.shape)

        flattened_idx = distance.argmin()
        distance.close()
        # NOTE: this has to correspond to NEMO nc file shape (y, x)
        ylat, xlon = np.unravel_index(flattened_idx, exp['toce'].shape)[2:]

        print('........ NEMO indices found:')
        print(xlon, ylat)

        xidx.append(xlon)
        yidx.append(ylat)
        print('....... checking list lengths:')
        print(len(xidx), len(yidx))

    # Select time_counter values closest to FLOAT timestamps
    # (done already at this point to get rid of everything that isn't necessary)
    tem = exp.toce.sel(time_counter = ftims, method='nearest', drop=True)
    tem = tem.drop_vars('time_centered')
    exp.close()

    # Create empty list to store NEMO vertical profiles to for each timestamp
    steps = []
    # List all NEMO timestamps
    tvals = tem.time_counter
    for ti, xi, yi in zip(tvals, xidx, yidx):
        # Select NEMO timestamp
        sel_step = tem.where(tem.time_counter==ti, drop=True)
        # Select corresponding lon and lat based on indices
        sel_step = sel_step.isel(x=xi, y=yi, drop=False)
        # Append to list
        steps.append(sel_step)

    # Form final dataframe from the list
    final = xr.concat(steps, dim='time_counter')
    print('final dataset')
    print(final)

    final.to_netcdf('float_with_nemo.nc', engine='netcdf4')
