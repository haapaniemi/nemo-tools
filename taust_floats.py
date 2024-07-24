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

def save_floats(list_float_ids, start_date, end_date):
    for float_no in floats:
        dat = argopy.DataFetcher().float(float_no)
        float_dat = dat.to_xarray()

        # Select time range
        fdata = float_dat.to_dataframe()
        time_mask = (fdata['TIME'] >= start_date) & (fdata['TIME'] <= end_date)
        fdata = fdata[time_mask]
        #print(fdata)
        #fdata = fdata.drop_duplicates('CYCLE_NUMBER')
        float_dat.to_netcdf('/scratch/project_2003984/Veera/run_nemo/results/taust_floats/float_'+ str(float_no) +'.nc')
    return

if __name__ == '__main__':

    # Define path to data
    # Select time range (currently used for FLOAT data only)

    print('starting...')
    floats = [6903711, 6903710, 6903707, 6903706, 6903703, 6903709, 6903708, 6903704]
    fetch_floats = False
    if fetch_floats:
        save_floats(floats, start_date, end_date)

    path2 = '/scratch/project_2003984/Veera/run_nemo/'
    #exp0 = xr.open_mfdataset(path2 + 'CTRL_OCT20/CTRL_1d_*_grid_T_*',\
    #                         engine='netcdf4', chunks='auto')
    exp0 = xr.open_mfdataset(path2 + 'TAUST_OCT20_Hs/TAUST_OCT20_Hs_1d_*_grid_T_*',\
                             engine='netcdf4', chunks='auto')

    ntims = exp0.time_counter.values
    ntim_min = pd.to_datetime(ntims[0]).strftime('%Y-%m-%d')
    ntim_max = pd.to_datetime(ntims[-1]).strftime('%Y-%m-%d')
    print('nemo times')
    print(ntim_min)
    print(ntim_max)

   # exp0 = xr.open_mfdataset(path2 + 'TAUST_OCT20_Hs/TAUST_OCT20_Hs_1d_*_grid_T_*',\
   #                          engine='netcdf4',chunks='auto')
   #  exp0 = xr.open_dataset(path2 +'CTRL_OCT20/CTRL_1d_20201002_20220101_grid_T_202010-202011.nc',engine='netcdf4')
    experiments = [['TAUST_OCT20_Hs_retry_2.0', exp0]]
    print('reading the exp data ok')

    for e in experiments:
        exp_name, exp = e[0], e[1]
        for float_no in floats:
            # If already in a file:
            fdata = xr.open_dataset(path2 + 'results/taust_floats/float_'+str(float_no)+'.nc')
            print('reading float data ok')

            # Drop duplicates --> only one time, lon, lat for each FLOAT CYCLE
            # (lat & lon same for each float cycle!)
            print(fdata)

            fdata = fdata.to_dataframe()
            fdata = fdata.reset_index()
            fdata = fdata.drop_duplicates('CYCLE_NUMBER')
            print('float data, dropped duplicates (dim = CYCLE_NUMBER):')

            fdata = fdata[['TIME', 'LONGITUDE', 'LATITUDE', 'TEMP', 'PRES']]
            print('float data reduced:')
            print(fdata)

            ftims_all = fdata['TIME'].values
            print('float timestamps')
            ftim_min = pd.to_datetime(ftims_all[0]).strftime('%Y-%m-%d')
            ftim_max = pd.to_datetime(ftims_all[-1]).strftime('%Y-%m-%d')
            print(ftim_min)
            print(ftim_max)

            start_times = [ntim_min, ftim_min]
            end_times = [ntim_max, ftim_max]

            # Define the overlapping timeframe
            start_date = max(start_times)
            end_date = min(end_times)
            print('selected timerange')
            print(start_date)
            print(end_date)

            # Select the same daterange from both nemo and float data
            exptim = exp.sel(time_counter=slice(start_date, end_date))
            time_mask = (fdata['TIME'] >= start_date) & (fdata['TIME'] <= end_date)
            fdata = fdata[time_mask]


            # Store float variables to data arrays for looping
            flats = fdata['LATITUDE'].values
            flons = fdata['LONGITUDE'].values
            ftims = fdata['TIME'].values
            print('float time stamps:')
            print(ftims)

            # Create empty lists to store the x and y indices tO
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

                distance = np.sqrt((exptim.nav_lon - lon)**2 + (exptim.nav_lat - lat)**2)
                print('........ distance set shape (only used for findind nearest point):')
                print(distance.shape)

                flattened_idx = distance.argmin()
                distance.close()
                # NOTE: this has to correspond to NEMO nc file shape (y, x)
                ylat, xlon = np.unravel_index(flattened_idx, exptim['toce'].shape)[2:]

                print('........ NEMO indices found:')
                print(xlon, ylat)

                xidx.append(xlon)
                yidx.append(ylat)
                print('....... checking list lengths:')
                print(len(xidx), len(yidx))

            # Select time_counter values closest to FLOAT timestamps
            # (done already at this point to get rid of everything that isn't necessary)
            tem = exptim.toce.sel(time_counter = ftims, method='nearest', drop=True)
            tem = tem.drop_vars('time_centered')
            exptim.close()
            #exp.close()

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

            final.to_netcdf('FLOAT_'+str(float_no)+'_'+exp_name +'_TEMP.nc', engine='netcdf4')
