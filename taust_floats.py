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

    path = '/scratch/project_2003984/Veera/run_nemo/'
   # exp0 = xr.open_mfdataset(path + 'CTRL_OCT20/CTRL_1d_*_grid_T_*',\
   #                          engine='netcdf4', chunks='auto')
   # exp0 = xr.open_mfdataset(path + 'TAUST_OCT20_Hs/TAUST_OCT20_Hs_1d_*_grid_T_*',\
   #                          engine='netcdf4', chunks='auto')
    exp0 = xr.open_mfdataset(path + 'TAUST-WA2_fit3_OCT20/TAUST-WA2_OCT20_1d_*_grid_T*',\
                             engine='netcdf4', chunks='auto')
    ntims = exp0.time_counter.values
    ntim_min = pd.to_datetime(ntims[0]).strftime('%Y-%m-%d')
    ntim_max = pd.to_datetime(ntims[-1]).strftime('%Y-%m-%d')
    print('nemo times')
    print(ntim_min)
    print(ntim_max)

    experiments = [['TAUST-WA2_fit3', exp0]]
    #experiments = [['TAUST_Hs',exp0]]
   # experiments = [['CTRL',exp0]]
    print('reading the exp data ok')

    var = 'sal'

    for e in experiments:
        exp_name, exp = e[0], e[1]
        for float_no in floats:
            # If already in a file:
            fdata = xr.open_dataset(path + 'results/taust_floats/float_'+str(float_no)+'.nc')
            print('reading float data ok')

            # Drop duplicates --> only one time, lon, lat for each FLOAT CYCLE
            # (lat & lon same for each float cycle!)
            print(fdata)

            fdata = fdata.to_dataframe()
            fdata = fdata.reset_index()
            fdata = fdata.drop_duplicates('CYCLE_NUMBER')
            print('float data, dropped duplicates (dim = CYCLE_NUMBER):')

            if var == 'temp':
                fdata = fdata[['TIME', 'LONGITUDE', 'LATITUDE', 'TEMP', 'PRES']]
            elif var == 'sal':
                fdata = fdata[['TIME', 'LONGITUDE', 'LATITUDE', 'PSAL', 'PRES']]

            print('float data reduced:')
            print(fdata)

            ftims_all = fdata['TIME'].values
            print('float timestamps')
            ftim_min = pd.to_datetime(ftims_all[0]).strftime('%Y-%m-%d')
            ftim_max = pd.to_datetime(ftims_all[-1]).strftime('%Y-%m-%d')

            start_times = [ntim_min, ftim_min]
            end_times = [ntim_max, ftim_max]

            # Define the overlapping timeframe
            start_date = max(start_times)
            end_date = min(end_times)

            # Select the same daterange from both nemo and float data
            exptim = exp.sel(time_counter=slice(start_date, end_date))
            time_mask = (fdata['TIME'] >= start_date) & (fdata['TIME'] <= end_date)
            fdata = fdata[time_mask]

            # Store float variables to data arrays for looping
            flats = fdata['LATITUDE'].values
            flons = fdata['LONGITUDE'].values
            ftims = fdata['TIME'].values
            print('float timestamps:', ftims)
            print('number of float timestamps', len(ftims))
            print('number of float lons', len(flons))
            print('number of float lats', len(flats))

            # Create empty lists to store the x and y indices tO
            xidx, yidx = [], []

            print('looping through float timestamps:')

            # Find the nearest indices for all of the time steps
            for i in range(len(ftims)):
                timestamp = ftims[i]
                print('........ current time stamp: ', timestamp)

                # Select float lon & lat for this timestamp
                lon = flons[i]
                lat = flats[i]
                print('........ float lon, lat: ', lon, lat)

                # Calculate the float's distance to nemo grid points
                distance = np.sqrt((exptim.nav_lon - lon)**2 + (exptim.nav_lat - lat)**2)
                print('........ distance set shape (only used for findind nearest point):')
                print(distance.shape)

                # Find the nearest grid point index
                flattened_idx = distance.argmin()
                distance.close()

                # NOTE: this has to correspond to NEMO nc file shape (y, x)
                if var == 'temp':
                    ylat, xlon = np.unravel_index(flattened_idx, exptim['toce'].shape)[2:]
                elif var == 'sal':
                    ylat, xlon = np.unravel_index(flattened_idx, exptim['soce'].shape)[2:]

                print('........ NEMO indices found:')
                print(xlon, ylat)

                xidx.append(xlon)
                yidx.append(ylat)
                print('....... checking list lengths:')
                print(len(xidx), len(yidx))

            print('final index list lengths')
            print(len(xidx), len(yidx))
            print('number of float timestamps')
            print(len(ftims))
            print(len(np.unique(ftims)))
            print(ftims)

            # Select time_counter values closest to FLOAT timestamps
            # (done already at this point to get rid of everything that isn't necessary)
            # These will then be assigned the float time stamps to avoid duplicates
            # (to solve issues with coarser temporal resolution)
            # --- should be fine up to 1d!
            if var == 'temp':
                tem = exptim.toce.sel(time_counter = ftims, method='nearest', drop=True)
            elif var == 'sal':
                tem = exptim.soce.sel(time_counter = ftims, method='nearest', drop=True)

            print('exp set after selecting timestamps')
            print(tem)
            tem = tem.drop_vars('time_centered')
            exptim.close()

            # Create empty list to store NEMO vertical profiles to for each timestamp
            steps = []
            print('tem old')
            print(tem.time_counter)
            tem = tem.assign_coords({'time_counter':ftims})
            #tem.time_counter.values = ftims
            print('tem new')
            print(tem.time_counter)
            tvals = tem.time_counter
            print('tem.time_counter values')
            print(tvals.values)
            print(tvals.shape)
            print(len(tvals))

            for ti, xi, yi in zip(tvals, xidx, yidx):
                sel_step = tem.where(tem.time_counter==ti, drop=True)
                sel_step = sel_step.isel(x=xi, y=yi, drop=True)

                print(sel_step.mean())
                print(f"Index: {i}, Time: {ti.values}, x: {xi}, y: {yi}")
                print(sel_step.time_counter.values)
                print(f"sel_step shape: {sel_step.shape}")
                steps.append(sel_step)

            print(f"Number of steps: {len(steps)}", len(steps))
            print(steps)

            # Form final dataframe from the list
            final = xr.concat(steps, dim='time_counter')
            print('final dataset')
            print(final)
            print(f"Final shape: {final.shape}")
            print(f'Number of float times: {len(ftims)}')

            if var == 'temp':
                final.to_netcdf(exp_name +'_' + str(float_no) +'_TEMP.nc', engine='netcdf4',\
                               encoding={'time_counter': {'dtype': 'i4'}})
            elif var == 'sal':
                final.to_netcdf(exp_name +'_' + str(float_no) + '_SAL.nc', engine='netcdf4',\
                               encoding={'time_counter': {'dtype': 'i4'}})

