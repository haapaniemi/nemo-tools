import xarray as xr

if __name__ == '__main__':
    # Read in model data
    path = '/scratch/project_2003984/Veera/run_nemo'

    # Read in 6h grid T data
#    exp1 = xr.open_mfdataset(path + '/TAUSTexp01.1/1h_files/TAUSTexp01.1_*grid_T*',\
#         chunks='auto',engine='netcdf4')
    exp2 = xr.open_mfdataset(path + '/TAUSTexp02.3/TAUSTexp02.3_1h_*grid_T*',\
         chunks='auto', engine='netcdf4')
    exp3 = xr.open_mfdataset(path + '/TAUSTexp02.4/TAUSTexp02.4_1h_*grid_T_*',\
         chunks='auto', engine='netcdf4')

    # [Names, files]
    list_files = [['TAUSTexp02.3_1h_', exp2],\
                  ['TAUSTexp02.4_1h_', exp3]] # ['TAUSTexp01.1_1h_', exp1],\
                  # ['TAUSTexp01.2_1h_', exp2],\


    # Specifying NEMO grid points to be selected
    nemo_points = [['Aarhus', 518, 460],\
                   ['Degerby', 883, 692],\
                   ['Forsmark', 806, 715],\
                   ['Frederikshavn', 530, 537],\
                   ['GoteborgTorshamnen', 574, 552],\
                   ['Hamina', 1128, 724],\
                   ['Hanko', 977, 680],\
                   ['Helsinki', 1048, 699],\
                   ['Hornbaek', 598, 457],\
                   ['Kemi', 1032, 1031],\
                   ['Kobenhavn', 604, 433],\
                   ['Kronstadt', 1221, 689],\
                   ['Kungsholmsfort', 710, 457],\
                   ['Kungsvik', 551, 630],\
                   ['LandsortNorra', 792, 617],\
                   ['Parnu', 1031, 593],\
                   ['Pietarsaari', 966, 913],\
                   ['Raahe', 1028, 970],\
                   ['Stockholm', 801, 651],\
                   ['Tallinn', 1041, 658],\
                   ['Viken', 602, 459],\
                   ['Visby', 807, 548],\
                   ['Warnemuende', 586, 342],\
                   ['Travemuende',542, 328]]

    # Loop through files in list_files
    for i in range(len(list_files)):
        # Select row
        sel_file = list_files[i]
        # Read name and NEMO data 
        site_name = sel_file[0]
        site_data = sel_file[1]

        # Loop though measurement points
        for j in range(len(nemo_points)):
            # Select next point from the list
            sel_point = nemo_points[j]
            # Select name and location of the site
            sel_point_name = sel_point[0]
            # Note: order is reversed compared to netcdf!
            sel_point_loc_x = sel_point[1]
            sel_point_loc_y = sel_point[2]

            selected_data = site_data.sel(x=sel_point_loc_x, y=sel_point_loc_y)
            selected_data.to_netcdf(site_name + sel_point_name +'.nc')
            print(site_name + sel_point_name + ' ready')

    print('Done')
