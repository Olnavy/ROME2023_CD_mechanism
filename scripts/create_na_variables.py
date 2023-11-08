'''
PhD 2 --- Create na variables for Figure 3 clusters
'''

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
# experiments = {'xoupa':"/nfs/see-fs-01_users/eeymr/database/xoupa", 
#                'tfgbi':"/nfs/see-fs-01_users/eeymr/database/tfgbi", 
#                'xoupk':"/nfs/see-fs-01_users/eeymr/database/xoupk",
#                'xoupf':"/nfs/see-fs-01_users/eeymr/database/xoupf",
#                'xouph':"/nfs/see-fs-01_users/eeymr/database/xouph", 
#                'xppbf':"/nfs/see-fs-01_users/eeymr/database/xppbf"}

experiments = {'xppbf':"/nfs/see-fs-01_users/eeymr/database/xppbf"}

ds_basin = xr.open_dataset("~/work/data/um/basin_hadcm3_glac1d_lgm.nc")

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

filt = util.ButterLowPass(order=1, fc=2*10**-3, fs=1, mult=2)


mld_zone = {}
stream_zone = {}
ocnt_zone, ocns_zone, ocnd_zone = {}, {}, {}

calc_mld = True
calc_stream = True
calc_ocntsd = False

for expt in experiments:

    if calc_mld:
        # Calculating MLD in NA zones
        print(f"__ Calculating MLD in NA zones for {expt}")

        temp_mld = xr.open_dataset(f"{experiments[expt]}/perso/{expt}.oceanmixedpf.annual.nc").mixLyrDpth_mm_uo
        
        if 'unspecified' in temp_mld.dims:
            temp_mld = temp_mld.isel(unspecified=0, drop=True)

        ds_zone = {}
        for zone in masks_na.keys():
            ds_zone[zone] = ((temp_mld*masks_na[zone]).mean(['longitude','latitude']).rename("mld")).expand_dims(zone=[zone])
            # ds_zone[zone].assign(mldf=filt.process((ds_zone[zone].mld).values))

        mld_zone[expt] = xr.concat(list(ds_zone.values()), dim='zone')


        # Saving the MLD in NA zones
        print(f"__ Saving MLD in NA zones for {expt}")

        sav_folder = f"{experiments[expt]}/perso/na_variables"
        util.makedir(sav_folder)
        
        util.nc_sav(mld_zone[expt], f"{sav_folder}/{expt}.mld.zone_na.annual")


    if calc_stream:
        # Calculating mean streamfunction in NA zones
        print(f"__ Calculating streamfunction in NA zones for {expt}")

        temp_stream = xr.open_dataset(f"{experiments[expt]}/perso/{expt}.streamFnpf01.annual.nc").streamFn_mm_uo
        
        if 'unspecified' in temp_stream.dims:
            temp_stream = temp_stream.isel(unspecified=0, drop=True)

        ds_zone = {}
        for zone in masks_na.keys():
            ds_zone[zone] = ((temp_stream*masks_na[zone]).mean(['longitude','latitude']).rename("stream")).expand_dims(zone=[zone])

        stream_zone[expt] = xr.concat(list(ds_zone.values()), dim='zone')


        # Saving the salinity clusters
        print(f"__ Saving streamfunction in NA zones for {expt}")

        sav_folder = f"{experiments[expt]}/perso/na_variables"
        util.makedir(sav_folder)
        
        util.nc_sav(stream_zone[expt], f"{sav_folder}/{expt}.stream.zone_na.annual")


    if calc_ocntsd:
        # Calculating mean ocean profiles in NA zones
        print(f"__ Calculating ocean profiles in NA zones for {expt}")

        print(f"____ Import Temperature")
        temp_ocnt = xr.open_dataset(f"{experiments[expt]}/time_series/{expt}.oceantemppg.annual.nc").temp_ym_dpth
        print(f"____ Import Salinity")
        temp_ocns = xr.open_dataset(f"{experiments[expt]}/time_series/{expt}.oceansalipg.annual.nc").salinity_ym_dpth*1000 + 35

        if os.path.exists(f"{experiments[expt]}/perso/{expt}.oceandenspg.annual.nc"):
            print(f"____ Import Density")
            temp_ocnd = xr.open_dataset(f"{experiments[expt]}/perso/{expt}.oceandenspg.annual.nc").density
        else:
            print(f"____ Calculate Density")
            temp_ocnd = util.density(temp_ocnt, temp_ocns).rename('density')
            util.nc_sav(temp_ocnd, f"{experiments[expt]}/perso/{expt}.oceandenspg.annual.nc")
    
        print(f"____ Calculate the na means")
        ds_zone_t, ds_zone_s, ds_zone_d = {}, {}, {}
        for zone in masks_na.keys():
            ds_zone_t[zone] = ((temp_ocnt*masks_na[zone]).mean(['longitude','latitude']).rename("temp_ym_dpth")).expand_dims(zone=[zone])
            ds_zone_s[zone] = ((temp_ocns*masks_na[zone]).mean(['longitude','latitude']).rename("salinity_ym_dpth")).expand_dims(zone=[zone])
            ds_zone_d[zone] = ((temp_ocnd*masks_na[zone]).mean(['longitude','latitude']).rename("density")).expand_dims(zone=[zone])

        ocnt_zone[expt] = xr.concat(list(ds_zone_t.values()), dim='zone')
        ocns_zone[expt] = xr.concat(list(ds_zone_s.values()), dim='zone')
        ocnd_zone[expt] = xr.concat(list(ds_zone_d.values()), dim='zone')


        # Saving the ocean profiles
        print(f"__ Saving ocean profiles in NA zones for {expt}")

        sav_folder = f"{experiments[expt]}/perso/na_variables"
        util.makedir(sav_folder)
        
        util.nc_sav(ocnt_zone[expt], f"{sav_folder}/{expt}.oceantemppg.zone_na.annual")
        util.nc_sav(ocns_zone[expt], f"{sav_folder}/{expt}.oceansalipg.zone_na.annual")
        util.nc_sav(ocnd_zone[expt], f"{sav_folder}/{expt}.oceandenspg.zone_na.annual")