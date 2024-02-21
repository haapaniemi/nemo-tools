"""
author: haapanie
"""
import xarray as xr
import matplotlib.pyplot as plt
import cmocean
import numpy as np
import pandas as pd
import gsw

def plot_model_comparison(ctrl, taust, var):
    if var == 'sal':
        cmap = cmocean.cm.haline #thermal
        vmin, vmax = 5, 6.6
        ctrl.soce.plot(x='time_counter', y='deptht',cmap=cmap, vmin=vmin, vmax=vmax)
    elif var == 'tem':
        cmap = cmocean.cm.thermal
        vmin, vmax = 4,18
        
    plt.title('CTRL')
    plt.ylim(150,0)
    plt.show()

    #TAUST
    taust.soce.plot(x='time_counter', y='deptht',cmap=cmap, vmin=vmin, vmax=vmax)
    plt.title('TAUST')
    plt.ylim(150,0)
    plt.show()

    # Difference between TAUST and CTRL
    dif = taust - ctrl
    dif.soce.plot(x='time_counter',y='deptht')
    plt.title('TAUST - CTRL')   # where values are < 0 --> CTRL > TAUST
    plt.ylim(150,0)             #                  > 0 --> CTRL < TAUST
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
    # Select only the important ones
    fdata = fdata[['TIME', 'PSAL','ABSSAL','PRES','DEPTH']]
    print(fdata)
    fdata = fdata.set_index(['TIME','DEPTH'])
    fdata = fdata.to_xarray()
    print(fdata)
    return fdata

if __name__ == '__main__':
    path = '/home/haapanie/'
    # -------------------------------------------------------------------------
    # Plotting model data only
    ctrl = xr.open_dataset(path+'FLOAT_6903710_SAL_CTRL_bfr_20210401-20210930.nc')
    taust = xr.open_dataset(path+'FLOAT_6903710_SAL_TAUSTexp01.2_20210401-20210930.nc')

    var = 'sal'
    plot_model_comparison(ctrl, taust, var)
    
    # -------------------------------------------------------------------------
    # Moving on with FLOAT data
    fdata = xr.open_dataset('float_6903710.nc')
    start_date = '2021-04-01'
    end_date = '2021-09-30'
    
    # Select time period
    fdata = sel_float_timerange(fdata, start_date, end_date)
    
    # Convert measurement data to model units
    fdata = convert_float_units(fdata)

    # Interpolate NEMO values to float depth points
    ctrl_new = ctrl.interp(deptht = fdata.DEPTH.values,method='linear')
    print('nemo data -----------------')
    print(ctrl_new)
    #print(ctrl_new.deptht.values)
    #print(fdata.PRES.values)
    #print(ctrl_new.deptht.values - fdata.PRES.values)

    plt.figure(figsize=(10,5))
    ctrl_new.soce.plot(x='time_counter',y='deptht')
    plt.ylim(150,0)
    plt.show()

    print(np.shape(ctrl_new.soce.values))
    print(np.shape(fdata.ABSSAL.values))

    dif = fdata.copy()
    print(dif)

    dif.ABSSAL.values = ctrl_new.soce.values - fdata.ABSSAL.values
    dif.ABSSAL.plot(x='TIME',y='DEPTH',vmin=-0.5, vmax=0.5,cmap=cmocean.cm.balance)
    plt.title('CTRL – ARGO')
    plt.ylim(150,0)
    plt.show()



    plt.figure(figsize=(8,5))
    x, y = np.meshgrid(dif['TIME'].values,dif['DEPTH'].values,indexing='ij')
    print(np.shape(x))
    print(np.shape(y))  
    plt.scatter(x,y, c=dif['ABSSAL'],cmap=cmocean.cm.balance,vmin=-1, vmax=1)   
    plt.colorbar()
    plt.ylim(150,0)
    plt.show()




    # -----------------------------------------------------------------------------
    taust = xr.open_dataset(path+'FLOAT_6903710_SAL_TAUST01_20210401-20210930.nc')  
    # Interpolate NEMO data to FLOAT depths
    # Calculate 
    
    taust_new = taust.interp(deptht=fdata.DEPTH.values)

    dif = fdata.copy()
    dif.ABSSAL.values = taust_new.soce.values - fdata.ABSSAL.values
    dif.ABSSAL.plot(x='TIME',y='DEPTH',vmin=-1, vmax=1,cmap=cmocean.cm.balance)
    plt.title('TAUST – ARGO')
    plt.ylim(150,0)
    plt.show()


#def plot_float_nemo_dif(data):
    plt.figure(figsize=(8,5))
    x, y = np.meshgrid(dif['TIME'].values, dif['DEPTH'].values, indexing='ij')
    print(np.shape(x))
    print(np.shape(y))
    plt.scatter(x,y, c=dif['ABSSAL'],cmap=cmocean.cm.balance,vmin=-1, vmax=1)
    plt.colorbar() 
    plt.ylim(150,0)
    plt.show()
