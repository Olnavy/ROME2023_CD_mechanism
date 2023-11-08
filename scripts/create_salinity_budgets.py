import xarray as xr
import numpy as np
import os 
from datetime import datetime
import cftime

import pylaeoclim_leeds.util_hadcm3 as util


'''
Compute the decadal salinity cubes and calculate their budgets for each simulation presented in PhD2
'''

print(f"\nRunning {os.path.basename(__file__)}")

# Global Variables
ds_volume = xr.open_dataset("/nfs/see-fs-01_users/eeymr/work/data/um/hadcm3.ocn_volume.nc")
n_lon_ocn = 13
experiments = ['xoupa', 'tfgbi', 'xoupk', 'xoupf', 'xouph', 'xppbf'] #PhD2


# Methods
def annual_to_decadal(expt):
    """
    Transform an annual time series to a decadal time series
    """

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    print(f"__ ({dt_string}) Annulal to decadal conversion for  {expt} - oceansalipg.annual")

    name_in = f"/nfs/see-fs-01_users/eeymr/database/{expt}/time_series/{expt}.oceansalipg.annual.nc"
    name_out = f"/nfs/see-fs-01_users/eeymr/database/{expt}/perso/{expt}.oceansalipg.decadal.nc"

    print(f"____ Opening {name_in}")
    ds_in = xr.open_dataset(name_in)

    print(f"____ Grouping by decade")
    ds_out = util.groupby_decade(
        ds_in.salinity_ym_dpth, shift_origin=False).rename("salinity_ym_dpth")
        
    print(f"____ Converting decade to cftime.Datetime360Day")
    ds_out = ds_out.rename(
        {'decade':'t'}).assign_coords(
            {'t':[cftime.Datetime360Day(decade*10,1,1) for decade in ds_out.decade]})

    util.nc_sav(ds_out, name_out)


def salinity_budget(expt):
    """
    Calculate the salinity budget of a salinity file
    """

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    print(f"__ ({dt_string}) Salinity budget for {expt} ")

    name_in = f"/nfs/see-fs-01_users/eeymr/database/{expt}/perso/{expt}.oceansalipg.decadal.nc"
    name_out = f"/nfs/see-fs-01_users/eeymr/database/{expt}/perso/salinity_budgets/{expt}.salinity_budget.decadal.nc"

    print(f"____ Opening {name_in}")
    ds_in = xr.open_dataset(name_in).salinity_ym_dpth*1000+35

    print(f"____ Calculating Budget")
    ds_out = (ds_in*ds_volume.volume).rename('salinity_budget')
    
    util.nc_sav(ds_out, name_out)

# ** Main method **

for expt in experiments:

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    
    print(f"\n({dt_string}) Processing {expt}\n---\n")
    
    # Create perso folder
    perso_folder =  f"/nfs/see-fs-01_users/eeymr/database/{expt}/perso/"
    util.makedir(perso_folder)

    # Annual to decade
    annual_to_decadal(expt)

    # Create budget folder
    budget_folder =  f"/nfs/see-fs-01_users/eeymr/database/{expt}/perso/salinity_budgets"
    util.makedir(budget_folder)

    # Calculate budget
    salinity_budget(expt)
