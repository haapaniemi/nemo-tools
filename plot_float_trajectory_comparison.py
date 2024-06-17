"""
author: haapanie
21/02/2024
"""
import xarray as xr
import matplotlib.pyplot as plt
import cmocean
import numpy as np
import pandas as pd
import gsw

smin, smax = 7,12
tmin, tmax = 4,18

ymax = 100

dsmin, dsmax = -2, 2

def plot_model_comparison(ctrl, taust, var):
    # Difference between TAUST and CTRL
    dif = taust - ctrl
    fig, axs = plt.subplots(1,3, sharex=True, sharey=True, figsize=(15,5))
    
    if var == 'sal':
        cmap = cmocean.cm.haline #thermal
        ctrl.soce.plot(x='time_counter', y='deptht',cmap=cmap,\
                       vmin=smin, vmax=smax, ax=axs[0])
        taust.soce.plot(x='time_counter', y='deptht',cmap=cmap,\
                         vmin=smin, vmax=smax, ax=axs[1])
        dif.soce.plot(x='time_counter',y='deptht',ax=axs[2])

    elif var == 'tem':
        cmap = cmocean.cm.thermal
        ctrl.toce.plot(x='time_counter', y='deptht',cmap=cmap,\
                       vmin=tmin, vmax=tmax, ax=axs[0])
        taust.toce.plot(x='time_counter', y='deptht',cmap=cmap,\
                         vmin=tmin, vmax=tmax, ax=axs[1])
        dif.toce.plot(x='time_counter',y='deptht',ax=axs[2])

            
    axs[0].set_ylim(150,0)
    axs[1].set_ylim(150,0)
    axs[2].set_ylim(150,0)

    axs[0].set_title('CTRL')
    axs[1].set_title('TAUST')
    axs[2].set_title('TAUST - CTRL')   # where values are < 0 --> CTRL > TAUST
                                        #                  > 0 --> CTRL < TAUST
    plt.show()
    return
def sel_float_timerange(fdata, start_date, end_date):
    fdata = fdata.to_dataframe()
    time_mask = (fdata['TIME'] >= start_date) & (fdata['TIME'] <= end_date)
    fdata = fdata[time_mask]
    return fdata
def convert_float_units(fdata):
    # From pressure to depth
    fdata['DEPTH'] = -gsw.conversions.z_from_p(fdata['PRES'].values, \
                                              np.mean(fdata['LATITUDE'].values))
    print('fdata----------------')
    print(fdata['DEPTH'])
    print(np.shape(fdata['DEPTH']))
    # From Practical Salinity to Absolute Salinity
    fdata['ABSSAL'] = gsw.conversions.SA_from_SP(fdata['PSAL'].values, \
                            fdata['DEPTH'].values, fdata['LONGITUDE'].values,\
                            fdata['LATITUDE'].values)
        
    fdata['CTEMP'] = gsw.conversions.CT_from_t(fdata['ABSSAL'], fdata['TEMP'],\
                                               fdata['PRES'])
        
    # Select only the important ones
    fdata = fdata[['TIME', 'PSAL','ABSSAL','PRES','DEPTH','CTEMP']]
    print(fdata)
    fdata = fdata.set_index(['TIME','DEPTH'])
    fdata = fdata.to_xarray()
    print(fdata)
    return fdata

def interp_NEMO_to_float(ctrl, taust, fdata, var):
    # Interpolate NEMO values to float depth points
    ctrl_new = ctrl.interp(deptht = fdata.DEPTH.values,method='linear')
    taust_new = taust.interp(deptht=fdata.DEPTH.values, method='linear')
    
    dif_ctrl = fdata.copy()
    dif_taust = fdata.copy()
      
    # Create mesh for plotting
    x, y = np.meshgrid(fdata['TIME'], fdata['DEPTH'],indexing='ij')
  
    # Sanity checking dimensions
    print(np.shape(x))
    print(np.shape(y)) 

    fig, axs = plt.subplots(1,5, sharex=True, sharey=True, figsize=(25,8))
    if var == 'sal':
        print(np.shape(ctrl_new.soce.values))
        print(np.shape(fdata.ABSSAL.values))
    
        smap = cmocean.cm.haline
        
        # Calculate differences
        dif_ctrl.ABSSAL.values = ctrl_new.soce.values - fdata.ABSSAL.values
        dif_taust.ABSSAL.values = taust_new.soce.values - fdata.ABSSAL.values
        
        # Plotting
        axs[0].scatter(x,y, c=fdata['ABSSAL'], cmap=smap,\
                       vmin=smin, vmax=smax)

        ctrl_new.soce.plot(x='time_counter',y='deptht', cmap=smap,\
                           vmin=smin, vmax=smax, ax=axs[1],add_colorbar=False)
        taust_new.soce.plot(x='time_counter',y='deptht',cmap=smap, \
                            vmin=smin, vmax=smax, ax=axs[2])
        
        axs[3].scatter(x,y, c=dif_ctrl['ABSSAL'], cmap=cmocean.cm.balance,\
                    vmin=dsmin, vmax=dsmax)
        a=axs[4].scatter(x,y, c=dif_taust['ABSSAL'], cmap=cmocean.cm.balance,\
                    vmin=dsmin, vmax=dsmax)
        plt.colorbar(a, ax=axs[4])
     #   dif_ctrl.ABSSAL.plot(x='TIME',y='DEPTH',vmin=smin, vmax=smax,\
     #                   cmap=cmocean.cm.balance,ax=axs[2])
     #   dif_taust.ABSSAL.plot(x='TIME', y='DEPTH', vmin=smin, vmax=smax,\
     #                         cmap=cmocean.cm.balance, ax=axs[3])
     
    elif var =='tem':
        tmin, tmax = -6,6
        tmap = cmocean.cm.thermal
        
        # Calculate differences
        dif_ctrl.CTEMP.values = ctrl_new.toce.values - fdata.CTEMP.values
        dif_taust.CTEMP.values = taust_new.toce.values - fdata.CTEMP.values
        
        # Plotting
        axs[0].scatter(x,y, c=fdata['CTEMP'], cmap=tmap)
        
        ctrl_new.toce.plot(x='time_counter',y='deptht', ax=axs[1],cmap=tmap,\
                           add_colorbar=False)
        taust_new.toce.plot(x='time_counter',y='deptht',ax=axs[2],cmap=tmap)
        
        axs[3].scatter(x,y, c=dif_ctrl['CTEMP'], cmap=cmocean.cm.balance,\
                    vmin=tmin, vmax=tmax)
        a=axs[4].scatter(x,y, c=dif_taust['CTEMP'], cmap=cmocean.cm.balance,\
                    vmin=tmin, vmax=tmax)
        plt.colorbar(a,ax=axs[4])
            
    axs[0].tick_params(labelrotation=45)
    axs[2].tick_params(labelrotation=45)
    axs[3].tick_params(labelrotation=45)
    axs[4].tick_params(labelrotation=45)

    axs[0].set_title('ARGO')
    axs[1].set_title('CTRL')
    axs[2].set_title('TAUST')
    axs[3].set_title('CTRL – ARGO')
    axs[4].set_title('TAUST – ARGO')

    axs[2].set_ylim(ymax,0)
    axs[3].set_ylim(ymax,0)    
    plt.show()
    
    return

if __name__ == '__main__':
    path = '/home/haapanie/'
    
    # FLOAT data
    fdata = xr.open_dataset('float_6903710.nc')
    start_date = '2021-04-01'
    end_date = '2021-09-30'
    
    # Select time period
    fdata = sel_float_timerange(fdata, start_date, end_date)
    # Convert measurement data to model units
    #   - using gsw package conversion functions
    #   - SA_from_SP
    #   - z_from_p
    #   - CT_from_t
    fdata = convert_float_units(fdata)
    
    # -------------------------------------------------------------------------
    # Plotting model data only
    ctrl = xr.open_dataset(path+'FLOAT_6903710_SAL_CTRL_bfr_20210401-20210930.nc')
    taust = xr.open_dataset(path+'FLOAT_6903710_SAL_TAUSTexp01.2_20210401-20210930.nc')
    # taust = xr.open_dataset(path+'FLOAT_6903710_SAL_TAUST01_20210401-20210930.nc')  
    
    var = 'sal'
    # Plot model data
    plot_model_comparison(ctrl, taust, var)
    
    # Plot NEMO & FLOAT comparison  
    interp_NEMO_to_float(ctrl, taust, fdata, var)
    

    # TEMP -----------
    ctrl = xr.open_dataset(path+'FLOAT_6903710_CTRL_bfr_20210401-20210930.nc')
    taust =  xr.open_dataset(path+'FLOAT_6903710_TAUSTexp01.2_20210401-20210930.nc')
    var='tem'
    plot_model_comparison(ctrl, taust, var)
    interp_NEMO_to_float(ctrl, taust, fdata, var)
    
    
    # Baltic Proper ARGO
    ctrl = xr.open_dataset(path + 'FLOAT_7900586_SAL_CTRL_bfr_20210401-20210930.nc')
    taust = xr.open_dataset(path + 'FLOAT_7900586_SAL_TAUST01.2_20210401-20210930.nc')
    fdata = xr.open_dataset('float_7900586.nc')
    
    fdata = sel_float_timerange(fdata, start_date, end_date)
    fdata = convert_float_units(fdata)
    
    var='sal'
    plot_model_comparison(ctrl, taust, var)
    interp_NEMO_to_float(ctrl, taust, fdata, var)

    var='tem'
    ctrl = xr.open_dataset(path + 'FLOAT_7900586_CTRL_bfr_20210401-20210930.nc')
    taust = xr.open_dataset(path + 'FLOAT_7900586_TAUSTexp01.2_20210401-20210930.nc')
    
    plot_model_comparison(ctrl, taust, var)
    interp_NEMO_to_float(ctrl, taust, fdata, var)
