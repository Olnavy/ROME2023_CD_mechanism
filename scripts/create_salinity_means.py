import xarray as xr
import numpy as np
import os 

import pylaeoclim_leeds.util_hadcm3 as util
from shapely.geometry import box, Point, Polygon


'''
Calculate the salinity clusters for each simulation presented in PhD2. Needs to be run after create_salinity_budgets.py
'''

print(f"Running {os.path.basename(__file__)}")

# Global variables
database = "/nfs/see-fs-01_users/eeymr/database"
ds_basin = xr.open_dataset("~/work/data/um/basin_hadcm3_glac1d_lgm.nc")
experiments = ['xoupa', 'tfgbi', 'xoupk', 'xoupf', 'xouph', 'xppbf']

depths = {'sfc':(5,667), 'int':(667,2000), 'deep':(2000,6000), 'tot':(None,None)}

masks = {}
masks['arc'] = ds_basin.arctic
masks['na'] = ds_basin.atlantic.sel(latitude=slice(44,79))
masks['tpa'] = ds_basin.atlantic.sel(latitude=slice(-33,44))
masks['sa'] = ds_basin.atlantic.sel(latitude=slice(-60,-33))
masks['pac'] = ds_basin.pacific
masks['so'] = ds_basin.southern
masks['idn'] = ds_basin.indian


masks_na={}
masks_na['gin'] = xr.concat([ds_basin.atlantic.sel(longitude=slice(340,360)).sel(latitude=slice(64,79)),
                          ds_basin.atlantic.sel(longitude=slice(0,20)).sel(latitude=slice(64,79))], 
                         dim='longitude')
masks_na['eur'] = xr.concat([ds_basin.atlantic.sel(longitude=slice(340,360)).sel(latitude=slice(44,64)),
                          ds_basin.atlantic.sel(longitude=slice(0,2)).sel(latitude=slice(44,64))], 
                         dim='longitude')
masks_na['irm'] = ds_basin.atlantic.sel(longitude=slice(316,339)).sel(latitude=slice(53,70))
masks_na['ls'] = ds_basin.atlantic.sel(longitude=slice(280,315)).sel(latitude=slice(53,79))
masks_na['spg'] = ds_basin.atlantic.sel(longitude=slice(280,339)).sel(latitude=slice(44,53))
masks_na['arc'] = ds_basin.arctic

# Chose between global and North Atlantic
mask_cluster = masks

# Import the salinity budgets
print('__ Importing ocean salinity')
salinity = {}
for expt in experiments: 
    salinity[expt] = xr.open_dataset(
        f"{database}/{expt}/perso/{expt}.oceansalipg.decadal.nc").salinity_ym_dpth*1000+35

salinity_means = {}

for expt in experiments:

    # Calculate salinity clusters
    print(f"__ Calculating salinity clusters for {expt}")
    salinity_means[expt] = {}
    ds_zone = {}

    for zone in mask_cluster.keys():
        print(f"__Processing {expt}-{zone}")
        ds_depth = {}
        
        for depth in ['sfc', 'int', 'deep','tot']:
            salinity_depth = salinity[expt].sel(depth_1=slice(*depths[depth]), drop=True)
            cluster = (salinity_depth*mask_cluster[zone]).weighted(
                salinity_depth.depth_1).mean(
                dim=['depth_1', 'longitude', 'latitude'])
            ds_depth[depth] = cluster.expand_dims(depth=[depth])

        ds_zone[zone] = xr.concat(list(ds_depth.values()), dim='depth').expand_dims(zone=[zone])

    salinity_means[expt] = (xr.concat(list(ds_zone.values()), dim='zone')).rename('salinity_means')


    # Saving the salinity clusters
    print(f"__ Saving salinity means for {expt}")

    sav_folder = f"/nfs/see-fs-01_users/eeymr/database/{expt}/perso/salinity_budgets"
    util.makedir(sav_folder)
    
    util.nc_sav(salinity_means[expt], f"{sav_folder}/{expt}.salinity_means.decadal")
