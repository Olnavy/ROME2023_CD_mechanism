import xarray as xr
import os 
import cftime

import pylaeoclim_leeds.util_hadcm3 as util


'''
Compute the decadal salinity cubes and calculate their budgets for each simulation presented in PhD2
'''

print(f"\nRunning {os.path.basename(__file__)}")


# Global Variables
depths = {'sfc':(5,667), 'int':(667,2000), 'deep':(2000,6000), 'tot':(None,None)}

ds_basin = xr.open_dataset("~/work/data/um/basin_hadcm3_glac1d_lgm.nc")

masks = {}

masks['arc'] = ds_basin.arctic
masks['na'] = ds_basin.atlantic.sel(latitude=slice(44,79))
masks['tpa'] = ds_basin.atlantic.sel(latitude=slice(-33,44))
masks['sa'] = ds_basin.atlantic.sel(latitude=slice(-60,-33))
masks['pac'] = ds_basin.pacific
masks['so'] = ds_basin.southern
masks['idn'] = ds_basin.indian


# Methods
budget_tot = lambda ds : (ds*masks[zone]).sum(
    dim=['latitude', 'longitude', 'depth_1'])*3600*24*360

budget_depth = lambda ds, depth_bounds : (ds*masks[zone]).sel(depth_1=slice(*depth_bounds)).sum(
    dim=['latitude', 'longitude', 'depth_1'])*3600*24*360

budget_sfc = lambda ds : (ds*masks[zone]).sum(
    dim=['latitude', 'longitude'])*3600*24*360


# Import advection volume dataset
print("__ Importing volume advection datasets")
data_folder = "/nfs/see-fs-01_users/eeymr/database/xoupk/perso/tendencies"

dsv = {
    'advx':xr.open_dataset(f"{data_folder}/xoupk.dsv_advx.decadal.nc").advx,
    'advy':xr.open_dataset(f"{data_folder}/xoupk.dsv_advy.decadal.nc").advy,
    'advz':xr.open_dataset(f"{data_folder}/xoupk.dsv_advz.decadal.nc").advz,
    'adv':xr.open_dataset(f"{data_folder}/xoupk.dsv_adv.decadal.nc").adv,
    'diffx':xr.open_dataset(f"{data_folder}/xoupk.dsv_diffx.decadal.nc").diffx,
    'diffy':xr.open_dataset(f"{data_folder}/xoupk.dsv_diffy.decadal.nc").diffy,
    'diffz':xr.open_dataset(f"{data_folder}/xoupk.dsv_diffz.decadal.nc").diffz,
    'conv':xr.open_dataset(f"{data_folder}/xoupk.dsv_conv.decadal.nc").conv,
    'fourier':xr.open_dataset(f"{data_folder}/xoupk.dsv_fourier.decadal.nc").fourier,
    'med':xr.open_dataset(f"{data_folder}/xoupk.dsv_med.decadal.nc").med,
    'gm':xr.open_dataset(f"{data_folder}/xoupk.dsv_gm.decadal.nc").gm,
    'ml':xr.open_dataset(f"{data_folder}/xoupk.dsv_ml.decadal.nc").ml,
    'robert':xr.open_dataset(f"{data_folder}/xoupk.dsv_robert.decadal.nc").robert,
    'sfc':xr.open_dataset(f"{data_folder}/xoupk.dsv_sfc.decadal.nc").sfc,
    'ice':xr.open_dataset(f"{data_folder}/xoupk.dsv_ice.decadal.nc").ice
}


# Calculate volume budget tendencies
print("__ Calculating volume budget tendencies")

ds_tdc = {}

for tdc in dsv.keys():
    ds_zone = {}

    for zone in masks.keys():
        print(f"____ Processing {tdc} - {zone}")
        ds_depth = {}

        if tdc not in ['sfc', 'ice']:
            for depth in ['sfc', 'int', 'deep']:
                ds_depth[depth] = (budget_depth(dsv[tdc], depth_bounds=depths[depth])).expand_dims(depth=[depth])
            ds_depth['tot'] = (budget_tot(dsv[tdc])).expand_dims(depth=['tot'])
        else:
            budget = budget_sfc(dsv[tdc])
            ds_depth['sfc'] = (budget).expand_dims(depth=['sfc'])
            ds_depth['int'] = (xr.zeros_like(budget)).expand_dims(depth=['int'])
            ds_depth['deep'] = (xr.zeros_like(budget)).expand_dims(depth=['deep'])
            ds_depth['tot'] = (budget).expand_dims(depth=['tot'])

            
        ds_zone[zone] = (xr.concat(list(ds_depth.values()), dim='depth')).expand_dims(zone=[zone])

    ds_tdc[tdc] = xr.concat(list(ds_zone.values()), dim='zone')

print(f"____ Calculating diff")
ds_tdc['diff'] = ds_tdc['diffx'] + ds_tdc['diffy'] + ds_tdc['diffz'] 

print(f"____ Calculating tot")
ds_tdc['tot'] = ds_tdc['adv'] + ds_tdc['diff'] + ds_tdc['conv'] + ds_tdc['fourier'] + ds_tdc['med'] + ds_tdc['gm'] + ds_tdc['ml'] + ds_tdc['robert'] + ds_tdc['sfc'] + ds_tdc['ice'] 

ds_out = xr.merge([ds_tdc[tdc].rename(tdc) for tdc in ds_tdc.keys()])    


# Converting cftime
print(f"____ Converting decade to cftime.Datetime360Day")
tdc_zone = ds_out.rename(
    {'decade':'t'}).assign_coords(
        {'t':[cftime.Datetime360Day(decade*10,1,1) for decade in ds_out.decade]})
    

# Saving tendencies zone
print("__ Saving volume budget tendencies")

util.nc_sav(tdc_zone, f"{data_folder}/xoupk.dsv.zone.decadal.nc")

