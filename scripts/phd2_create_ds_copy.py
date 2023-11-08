import xarray as xr
import os
import pylaeoclim_leeds.util_hadcm3 as util
from datetime import datetime
import cftime


print (f"\n--------\nRunning {os.path.basename(__file__)}\n--------\n")


# Global Variables
month_name={1:'jn', 2:'fb', 3:'mr', 4:'ap', 5:'my', 6:'jn', 7:'jl', 8:'ag', 9:'sp', 10:'ot', 11:'nv', 12:'dc'}
ds_volume = xr.open_dataset("/nfs/see-fs-01_users/eeymr/work/data/um/hadcm3.ocn_volume.nc")
# For non-extended dataset, set n_lon_ocn to 0
n_lon_ocn = 0


# Experiments list 
# experiments = ['xpfje', 
#     'xppbe', 'xppbf']
experiments = ['xpfje', 'xppbf']


# Methods

def create_folder(expt):
    """
    Create perso time series folder if not already existing
    """
    if not os.path.exists(name_out):
        print(f"__ Creating {expt}/perso/")
        os.mkdir(f"/nfs/see-fs-01_users/eeymr/database/{expt}/perso/")


def annual_to_decadal(expt, variable):
    """
    Transform an annual time series to a decadal time series
    """

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    print(f"__ ({dt_string}) Annulal to decadal conversion for {expt} - {variable}")

    name_in = f"/nfs/see-fs-01_users/eeymr/database/{expt}/{variable['folder']}/{expt}.{variable['name']}.nc"
    name_out = name_in.replace('annual','decadal').replace(variable['folder'],'perso')

    print(f"____ Opening {name_in}")
    ds_in = xr.open_dataset(name_in)

    print(f"____ Grouping by decade")
    if 'longitude_extended' not in ds_in.attrs or ds_in.longitude_extended == 'False':
        ds_out = util.groupby_decade(
            util.extend_lon(
                ds_in[variable['var_name']], n_lon_ocn, lon_name=variable['lon']), 
                shift_origin=False)
        ds_out.attrs['longitude_extended'] = 'True' if n_lon_ocn!=0 else 'False'
    else:
        ds_out = util.groupby_decade(ds_in[variable['var_name']], shift_origin=False)
        
    print(f"____ Converting decade to cftime.Datetime360Day")
    ds_out = ds_out.rename(
        {'decade':'t'}).assign_coords(
            {'t':[cftime.Datetime360Day(decade*10,1,1) for decade in ds_out.decade]})

    print(f"____ Saving at {name_out}")
    if os.path.exists(name_out):
        os.remove(name_out)
    ds_out.to_netcdf(name_out)


def monthly_to_annual(expt, variable):
    """
    Transform a monthly time series into an annual time series
    """

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    print(f"__ ({dt_string}) Annulal to decadal conversion for {expt} - {variable['name']}")

    name_in = f"/nfs/see-fs-01_users/eeymr/database/{expt}/{variable['folder']}/{expt}.{variable['name']}.nc"
    name_out = name_in.replace('monthly','annual').replace(variable['folder'],'perso')

    print(f"____ Opening {name_in}")
    ds_in = xr.open_dataset(name_in)

    print(f"____ Grouping by decade")
    if 'longitude_extended' not in ds_in.attrs or ds_in.longitude_extended == 'False':
        ds_out = util.extend_lon(
            ds_in[variable['var_name']].groupby('t.year').mean('t'), 
            n_lon_ocn, lon_name=variable['lon'])
        ds_out.attrs['longitude_extended'] = 'True' if n_lon_ocn!=0 else 'False'
    else:
        ds_out = ds_in[variable['var_name']].groupby('t.year').mean('t')
    
    print(f"____ Converting year to cftime.Datetime360Day")
    ds_out = ds_out.rename(
        {'year':'t'}).assign_coords(
            {'t':[cftime.Datetime360Day(year,1,1) for year in ds_out.year]})

    print(f"____ Saving at {name_out}")
    if os.path.exists(name_out):
        os.remove(name_out)
    ds_out.to_netcdf(name_out)


def monthly_to_month(expt, variable, month):
    """
    Transform a monthly time series into a month time series
    month: integer code for month
    """

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    print(f"__ ({dt_string}) Annulal to decadal conversion for {expt} - {variable['name']}")

    name_in = f"/nfs/see-fs-01_users/eeymr/database/{expt}/{variable['folder']}/{expt}.{variable['name']}.nc"
    name_out = name_in.replace('monthly',month_name[month]).replace(variable['folder'],'perso')

    print(f"____ Opening {name_in}")
    ds_in = xr.open_dataset(name_in)

    print(f"____ Grouping by {month_name[month]}")
    if 'longitude_extended' not in ds_in.attrs or ds_in.longitude_extended == 'False':
        ds_out = util.extend_lon(
            ds_in[variable['var_name']].sel(t=ds_in.t.dt.month.isin([month])), 
            n_lon_ocn, lon_name=variable['lon'])
        ds_out.attrs['longitude_extended'] = 'True' if n_lon_ocn!=0 else 'False'

    else:
        ds_in[variable['var_name']].sel(t=ds_in.t.dt.month.isin([month]))

    print(f"____ Saving at {name_out}")
    if os.path.exists(name_out):
        os.remove(name_out)
    ds_out.to_netcdf(name_out)


def salinity_budget(expt, variable):
    """
    Calculate the salinity budget of a salinity file
    """

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    print(f"__ ({dt_string}) Salinity budget for {expt} - {variable['name']}")

    name_in = f"/nfs/see-fs-01_users/eeymr/database/{expt}/{variable['folder']}/{expt}.{variable['name']}.nc"
    name_out = name_in.replace('oceansalipg','salinity_budget').replace(variable['folder'],'perso')

    print(f"____ Opening {name_in}")
    ds_in = xr.open_dataset(name_in)[variable['var_name']]*1000+35

    print(f"____ Calculating Budget")
    if 'longitude_extended' in ds_in.attrs and ds_in.longitude_extended == 'True':
        ds_out = (ds_in*util.extend_lon(
            ds_volume.volume, n_lon_ocn, lon_name='longitude')).rename('budget')
        ds_out.attrs['longitude_extended'] = 'True'
    else:
        ds_out = (ds_in*ds_volume.volume).rename('budget')
        ds_out.attrs['longitude_extended'] = 'False'

    
    print(f"____ Saving at {name_out}")
    if os.path.exists(name_out):
        os.remove(name_out)
    ds_out.to_netcdf(name_out)


def density_calculation(expt, variables):
    """
    Calculate the density from a temperature and salinity file
    """

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    print(f"__ ({dt_string}) Calculate density for {expt} - {variables['T']['name']} and {variables['S']['name']}")

    temperature_in = f"/nfs/see-fs-01_users/eeymr/database/{expt}/{variable['T']['folder']}/{expt}.{variable['T']['name']}.nc"
    salinity_in = f"/nfs/see-fs-01_users/eeymr/database/{expt}/{variable['S']['folder']}/{expt}.{variable['S']['name']}.nc"
    name_out = salinity_in.replace('oceansalipg','oceandenspg').replace(variable['S']['folder'],'perso')

    print(f"____ Opening {salinity_in} and {temperature_in}")
    ds_temperature_in = xr.open_dataset(temperature_in)[variable['T']['var_name']]
    ds_salinity_in = xr.open_dataset(salinity_in)[variable['S']['var_name']]*1000+35

    print(f"____ Calculating density")
    ds_out = util.density(ds_temperature_in, ds_salinity_in).rename('density')
    
    print(f"____ Saving at {name_out}")
    if os.path.exists(name_out):
        os.remove(name_out)
    ds_out.to_netcdf(name_out)


# ** Main method **

for expt in experiments:

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    
    print(f"({dt_string}) Processing {expt} ")
    
    # Monthly to annual
    # monthly_to_annual_in = [
    #     {'lon':'longitude','var_name':'mixLyrDpth_mm_uo', 'folder':'time_series', 'name':'oceanmixedpf.monthly'},
    #     {'lon':'longitude','var_name':'streamFn_mm_uo', 'folder':'time_series', 'name':'streamFnpf01.monthly'},
    #     {'lon':'longitude','var_name':'p_mm_msl', 'folder':'time_series', 'name':'mslp.monthly'}
    #     ]


    # for variable in monthly_to_annual_in:
    #     monthly_to_annual(expt, variable)

    # Monthly to month
    # monthly_to_march_in = [
    #     {'lon':'longitude','var_name':'iceconc_mm_srf', 'folder':'time_series', 'name':'iceconc.monthly'}
    #     ]
    # for variable in monthly_to_march_in:
    #     monthly_to_month(expt, variable, 3)

    # monthly_to_september_in = [
    #     {'lon':'longitude','var_name':'iceconc_mm_srf', 'folder':'time_series', 'name':'iceconc.monthly'}
    #     ]
    # for variable in monthly_to_september_in:
    #     monthly_to_month(expt, variable, 9)

    # Annual to decade
    annual_to_decade_in = [
        {'lon':'longitude','var_name':'salinity_ym_dpth', 'folder':'time_series', 'name':'oceansalipg.annual'},
        {'lon':'longitude','var_name':'temp_ym_dpth', 'folder':'time_series', 'name':'oceantemppg.annual'},
        ]
    #     {'lon':'longitude_1','var_name':'ucurrTot_ym_dpth', 'folder':'time_series', 'name':'oceanuvelpg.annual'},
    #     {'lon':'longitude_1','var_name':'vcurrTot_ym_dpth', 'folder':'time_series', 'name':'oceanvvelpg.annual'}, 
    #     {'lon':'longitude','var_name':'p_mm_msl', 'folder':'perso', 'name':'mslp.annual'}
    #     ]

    for variable in annual_to_decade_in:
        annual_to_decadal(expt, variable)


    # Calculate budget
    salinity_budget_in = [
        {'lon':'longitude','var_name':'salinity_ym_dpth', 'folder':'perso', 'name':'oceansalipg.decadal'},
        ]

    for variable in salinity_budget_in:
        salinity_budget(expt, variable)


    # Calculate density
    density_in = [{
        'T':{'lon':'longitude','var_name':'temp_ym_dpth', 'folder':'perso', 'name':'oceantemppg.decadal'},
        'S':{'lon':'longitude','var_name':'salinity_ym_dpth', 'folder':'perso', 'name':'oceansalipg.decadal'}
        }
        ]

    for variable in density_in:
        density_calculation(expt, variable)




