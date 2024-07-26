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

#smin, smax = 7,12
#tmin, tmax = 4,18

smin, smax = 5,15
tmin, tmax = 4,18

ymax = 250

dsmin, dsmax = -3, 3

ds_dif_min, ds_dif_max = -0.5, 0.5
dt_dif_min, dt_dif_max = -3, 3

def plot_model_comparison(ctrl, taust, var, fdata, title1, title2, float_id):
    print('here -----------------------')
    if float_id in [6903707,6903710,6903711]:
        smin, smax = 2, 7
        print('limits checked')
    else:
        smin, smax = 5,15
        
    if float_id in [6903707]:
        title = 'Bothnian Bay'
    elif float_id in [6903710, 6903711]:
        title = 'Bothnian Sea'
    
    # Difference between TAUST and CTRL
    dif = taust - ctrl
    fig, axs = plt.subplots(1,4, sharex=True, sharey=True, figsize=(17,5))
    
    ftemp = fdata.CTEMP.values
    fsal = fdata.ABSSAL.values
    # Create mesh for plotting
    x, y = np.meshgrid(fdata['TIME'], fdata['DEPTH'],indexing='ij')
    
    if var == 'sal':
        cmap = cmocean.cm.haline #thermal
        ctrl.soce.plot(x='time_counter', y='deptht',cmap=cmap,\
                       vmin=smin, vmax=smax, ax=axs[0], add_colorbar=False)
        a = taust.soce.plot(x='time_counter', y='deptht',cmap=cmap,\
                         vmin=smin, vmax=smax, ax=axs[1], add_colorbar=False)
        dif.soce.plot(x='time_counter',y='deptht',ax=axs[2],\
                      vmin=ds_dif_min, vmax=ds_dif_max, cmap=cmocean.cm.balance)
        
        axs[3].scatter(x,y, c=fsal, cmap=cmap, vmin=smin, \
                      vmax = smax)
        plt.colorbar(a, ax=axs[1])

    elif var == 'tem':
        cmap = cmocean.cm.thermal
        ctrl.toce.plot(x='time_counter', y='deptht',cmap=cmap,\
                       vmin=tmin, vmax=tmax, ax=axs[0], add_colorbar=False)
        a = taust.toce.plot(x='time_counter', y='deptht',cmap=cmap,\
                         vmin=tmin, vmax=tmax, ax=axs[1], add_colorbar=False)
        dif.toce.plot(x='time_counter',y='deptht',ax=axs[2],cmap=\
                      cmocean.cm.balance, vmin=dt_dif_min, vmax=dt_dif_max)
        axs[3].scatter(x,y, c=ftemp, cmap=cmap, vmin=tmin, vmax=tmax)
            
        plt.colorbar(a, ax=axs[1])

    plt.suptitle(float_id)
    axs[0].set_ylim(ymax,0)
    axs[1].set_ylim(ymax,0)
    axs[2].set_ylim(ymax,0)
    axs[3].set_ylim(ymax,0)
    
    axs[0].tick_params(labelrotation=45)
    axs[1].tick_params(labelrotation=45)
    axs[2].tick_params(labelrotation=45)
    axs[3].tick_params(labelrotation=45)

    axs[0].set_xlabel(' ')
    axs[1].set_xlabel(' ')
    axs[2].set_xlabel(' ')
    axs[3].set_xlabel(' ')
   
    axs[0].set_ylabel(' ')
    axs[0].set_ylabel(' ')
    axs[0].set_ylabel(' ')
    axs[0].set_ylabel(' ')


    axs[0].set_title('A: '+ title1)
    axs[1].set_title('B: '+ title2)
    axs[2].set_title('B - A')   # where values are < 0 --> CTRL > TAUST
                                        #                  > 0 --> CTRL < TAUST
    plt.show()
    return
def sel_float_timerange(fdata, start_date, end_date):
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

def interp_NEMO_to_float(ctrl, taust, fdata, var, title1, title2, float_id):
    
    if float_id in ['6903707','6903710','6903711']:
        smin, smax = 2, 7
    else:
        smin, smax = 5,15

    # -------------------------------------------------------------------------
    print('in interp_NEMO_to_float')
    print(fdata)
    fdata = fdata.to_dataframe().reset_index()
    print(fdata)
#    print(fdata.index.get_level_values[0])
    fdata = fdata.groupby([fdata.TIME.dt.date, fdata.DEPTH]).mean()
    fdata = fdata.to_xarray()
    
    
    print(ctrl)
    ctrl = ctrl.to_dataframe().reset_index()
    print(ctrl)
    ctrl = ctrl.groupby([ctrl.time_counter.dt.date, ctrl.deptht]).mean()
    print(ctrl)
    ctrl = ctrl.to_xarray()
    
    taust = taust.to_dataframe().reset_index()
    print(taust)
    taust = taust.groupby([taust.time_counter.dt.date, taust.deptht]).mean()
    print(taust)
    taust = taust.to_xarray()
    
    ftims_all = fdata['TIME'].values
    ftim_min = pd.to_datetime(ftims_all[0])#.strftime('%Y-%m-%d')
    ftim_max = pd.to_datetime(ftims_all[-1])#.strftime('%Y-%m-%d')
    
    print(ftim_min)
    print(ftim_max)
    
    ctrl = ctrl.sel(time_counter=slice(ftim_min, ftim_max))
    taust = taust.sel(time_counter=slice(ftim_min, ftim_max))
    
    # -----------------------------------------------------------------------
    
    
    # Interpolate NEMO values to float depth points
    ctrl_new = ctrl.interp(deptht = fdata.DEPTH.values,method='linear')
    taust_new = taust.interp(deptht=fdata.DEPTH.values, method='linear')
    
    dif_ctrl = fdata.copy()
    dif_taust = fdata.copy()
      
    # Create mesh for plotting
    x, y = np.meshgrid(fdata['TIME'], fdata['DEPTH'],indexing='ij')
  
    # Sanity checking dimensions
    print('float data dimensions mesh')
    print('x', np.shape(x))
    print('y', np.shape(y)) 

    fig, axs = plt.subplots(1,5, sharex=True, sharey=True, figsize=(25,8))
    plt.suptitle(float_id)
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
        tmin, tmax = -4,4
        tmap = cmocean.cm.thermal
   
        print('dimensions before calculating differences')
        print('nemo 1', np.shape(ctrl_new.toce.values))
        print('nemo 2', np.shape(taust_new.toce.values))
        print('float', np.shape(fdata.CTEMP.values))
    
       # smap = cmocean.cm.haline
        
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
            
    plt.title(float_id)
    axs[1].set_ylabel(' ')
    axs[2].set_ylabel(' ')
    axs[1].set_xlabel(' ')
    axs[2].set_xlabel(' ')
    
    axs[0].tick_params(labelrotation=45)
    axs[1].tick_params(labelrotation=45)
    axs[2].tick_params(labelrotation=45)
    axs[3].tick_params(labelrotation=45)
    axs[4].tick_params(labelrotation=45)
    
    print('here-------------------')
    dif_vals_ctrl = dif_ctrl['CTEMP'].values
    dif_vals_ctrl = dif_vals_ctrl[~np.isnan(dif_vals_ctrl)]
    print(dif_vals_ctrl)

#    print(dif_ctrl)
#    print(dif_taust)
    
    dif_vals_taust = dif_taust['CTEMP'].values
    dif_vals_taust = dif_vals_taust[~np.isnan(dif_vals_taust)]
    print(dif_vals_taust)
    
    
    rmse_1 = np.sqrt((dif_vals_ctrl)**2).mean()
    rmse_2 = np.sqrt((dif_vals_taust)**2).mean()
    
    print(rmse_1)
    print(rmse_2)
    

    axs[0].set_title('ARGO')
    axs[1].set_title('A: '+ title1)
    axs[2].set_title('B: '+ title2)
    axs[3].set_title('A – ARGO: RMSE = {:.2f}'.format(rmse_1))
    axs[4].set_title('B – ARGO: RMSE = {:.2f}'.format(rmse_2))

    axs[2].set_ylim(ymax,0)
    axs[3].set_ylim(ymax,0)    
    plt.show()
    
    return
def define_common_times(fdata, exp0, exp1, exp2):
    
  #  exp0 = exp0.drop_duplicates('time_counter')
  #  exp1 = exp1.drop_duplicates('time_counter')
  #  exp2 = exp2.drop_duplicates('time_counter')
   # fdata = fdata.drop_duplicates('TIME')
    
    n0tims = exp0.time_counter.values
    n0tim_min = pd.to_datetime(n0tims[0])#.strftime('%Y-%m-%d')
    n0tim_max = pd.to_datetime(n0tims[-1])#.strftime('%Y-%m-%d')
    
    n1tims = exp1.time_counter.values
    n1tim_min = pd.to_datetime(n1tims[0])#.strftime('%Y-%m-%d')
    n1tim_max = pd.to_datetime(n1tims[-1])#.strftime('%Y-%m-%d')
    
    n2tims = exp2.time_counter.values
    n2tim_min = pd.to_datetime(n2tims[0])#.strftime('%Y-%m-%d')
    n2tim_max = pd.to_datetime(n2tims[-1])#.strftime('%Y-%m-%d')
    
    ftims_all = fdata['TIME'].values
  #  print('float timestamps')
    ftim_min = pd.to_datetime(ftims_all[0])#.strftime('%Y-%m-%d')
    ftim_max = pd.to_datetime(ftims_all[-1])#.strftime('%Y-%m-%d')
 #   print(ftim_min)
  #  print(ftim_max)

#    start_times = [n0tim_min, n1tim_min, n2tim_min, ftim_min]
#    end_times = [n0tim_max, n1tim_max, n2tim_max, ftim_max]
    
    
    start_times = [n0tim_min, n2tim_min, ftim_min]
    end_times = [n0tim_max, n2tim_max, ftim_max]
    
    # Define the overlapping timeframe
    start_date = max(start_times)
    end_date = min(end_times)
    
    print('selected time range')
    print(start_date)
    print(end_date)
    
    fdata = sel_float_timerange(fdata, start_date, end_date)
    exp0 = exp0.sel(time_counter=slice(start_date, end_date))
    exp1 = exp1.sel(time_counter=slice(start_date, end_date))
    exp2 = exp2.sel(time_counter=slice(start_date, end_date))
    
    print(np.unique(fdata.TIME.values))
    print(np.unique(exp0.time_counter.values))
    
    return fdata, exp0, exp1, exp2

if __name__ == '__main__':
    path = '/home/haapanie/taust_manuscript/argo_plots/'
   
    float_ids = [6903711, 6903710, 6903707, 6903706, 6903703, 6903709, 6903708, 6903704]
    var = 'temp'

  #  float_ids = [6903711
  #  float_ids = [6903708]
    
    for float_id in float_ids:
        
        if var == 'sal':
            varname = 'SAL'
        elif var == 'temp':
            varname = 'TEMP'
        
        print('currently considering: ', float_id)
        # Model data
        exp0 = xr.open_dataset(path+'CTRL_'+str(float_id)+'_'+varname+'.nc')
        exp1 = xr.open_dataset(path+'TAUST_Hs_'+str(float_id)+'_'+varname+'.nc')
        exp2 = xr.open_dataset(path+'TAUST-WA2_fit3_'+str(float_id)+'_'+varname+'.nc')
        # FLOAT data
        fdata = xr.open_dataset(path+'float_'+str(float_id)+'.nc')
        fdata = fdata.to_dataframe()
       
        
        fdata, exp0, exp1, exp2 = define_common_times(fdata, exp0, exp1, exp2)

        print('---- after define_common_times')
        print(fdata)
        print(exp0)
        print(exp1)
        print(exp2)

        print('lengths')
        print('float times', len(np.unique(fdata.TIME.values)))
        print('exp0 times', len(np.unique(exp0.time_counter.values)))
        print('exp1 times', len(np.unique(exp1.time_counter.values)))
        print('exp2 times', len(np.unique(exp2.time_counter.values))) 
        # Select time period
        #fdata = sel_float_timerange(fdata, start_date, end_date)
        
        # Convert measurement data to model units
        #   - using gsw package conversion functions
        #   - SA_from_SP
        #   - z_from_p
        #   - CT_from_t
        fdata = convert_float_units(fdata)
#        print(fdata)
#        print(exp0)
#        print(exp1)
#        print(exp2)
        
#        print(k)
        # taust = xr.open_dataset(path+'FLOAT_6903710_SAL_TAUST01_20210401-20210930.nc')  
 #       ctrl = ctrl.sel(time_counter=slice('2021-04-01', '2021-08-31'))
 #       taust = ctrl.sel(time_counter=slice('2021-04-01', '2021-08-31'))

        plotting = True
        if plotting:
            # Plot model data
       #     plot_model_comparison(exp0, exp1, var, fdata, 'CTRL', 'TAUST_1',float_id)
            plot_model_comparison(exp0, exp2, var, fdata, 'CTRL', 'TAUST_2',float_id)
       #     plot_model_comparison(exp1, exp2, var, fdata, 'TAUST_1', 'TAUST_2',float_id)
            
            
      #  interp_NEMO_to_float(exp0, exp1, fdata, var, 'CTRL', 'TAUST_1', float_id)
        interp_NEMO_to_float(exp0, exp2, fdata, var, 'CTRL', 'TAUST_2', float_id)
     #   interp_NEMO_to_float(exp1, exp2, fdata, var, 'TAUST_1', 'TAUST_2', float_id)


       # print(tuple(exp0.dims[d] for d in ['time_counter', 'deptht']))
       ## print(fdata.shape)
        #print(fdata)
        #print(exp0)
        #exp0_new = exp0.interp(deptht = fdata.DEPTH.values,method='linear')
        #exp1_new = exp1.interp(deptht = fdata.DEPTH.values,method='linear')

        #print(exp0_new)
        
        #plot_model_comparison(exp0_new, exp1_new, var, fdata, 'CTRL', 'TAUST_1',float_id)

        #taust_new = taust.interp(deptht=fdata.DEPTH.values, method='linear')
            
        # Plot NEMO & FLOAT comparison      #

        # TEMP -----------
    #    ctrl = xr.open_dataset(path+'FLOAT_6903710_CTRL_bfr_20210401-20210930.nc')
     #   taust =  xr.open_dataset(path+'FLOAT_6903710_TAUSTexp01.2_20210401-20210930.nc')
     #   var='tem'
        
      #  plot_model_comparison(ctrl, taust, var, 'CTRL_bfr', 'TAUSTexp01.2')
     #   interp_NEMO_to_float(ctrl, taust, fdata, var, 'CTRL_bfr', 'TAUSTexp01.2')
    
