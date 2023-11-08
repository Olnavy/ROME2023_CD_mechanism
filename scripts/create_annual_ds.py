'''
PhD 2 --- Create annual datasets
'''

import xarray as xr
import numpy as np
import os 
from datetime import datetime

import pylaeoclim_leeds.util_hadcm3 as util
from shapely.geometry import box, Point, Polygon
import cftime

'''
Compute annual datasets from monthly datasets
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

variables = {'oceanmixedpf':'mixLyrDpth_mm_uo',
             'streamFnpf01':'streamFn_mm_uo'}


calc_mld = True
calc_stream = True


def monthly_to_annual(expt, var):
    """
    Transform a monthly time series into an annual time series
    """

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    print(f"__ ({dt_string}) Annulal to decadal conversion for {expt} - {var}")

    name_in = f"{experiments[expt]}/time_series/{expt}.{var}.monthly.nc"

    print(f"____ Opening {name_in}")
    ds_in = xr.open_dataset(name_in)

    print(f"____ Grouping by year")
    ds_out = ds_in[variables[var]].groupby('t.year').mean('t')
    
    print(f"____ Converting year to cftime.Datetime360Day")
    ds_out = ds_out.rename(
        {'year':'t'}).assign_coords(
            {'t':[cftime.Datetime360Day(year,1,1) for year in ds_out.year]})

    return ds_out

for expt in experiments:

    if calc_mld:
        # Calculating MLD
        print(f"__ Calculating MLD annual means for {expt}")

        mld = monthly_to_annual(expt, 'oceanmixedpf')

        # Saving the MLD
        print(f"__ Saving MLD annual means for {expt}")

        sav_folder = f"{experiments[expt]}/perso/"
        util.makedir(sav_folder)
        
        util.nc_sav(mld.isel(unspecified=0, drop=True), f"{sav_folder}/{expt}.oceanmixedpf.annual")


    if calc_stream:
        # Calculating streamfunction
        print(f"__ Calculating streamfunction annual means for {expt}")


        stream = monthly_to_annual(expt, 'streamFnpf01')

        # Saving the streamfunction
        print(f"__ Saving streamfunction annual means for {expt}")

        sav_folder = f"{experiments[expt]}/perso/"
        util.makedir(sav_folder)
        
        util.nc_sav(stream.isel(unspecified=0, drop=True), f"{sav_folder}/{expt}.streamFnpf01.annual")

