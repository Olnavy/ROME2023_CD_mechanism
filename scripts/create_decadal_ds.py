import xarray as xr
import pylaeoclim_leeds.util_hadcm3 as util
import os

expt = 'xoupk'
dataset_folder = "/nfs/see-fs-01_users/eeymr/database/xoupk"
n_lon = 13

na_lon = (260,380)

variables = {'S':{'lon':'longitude','var_name':'salinity_ym_dpth', 'folder':'time_series', 'name':'xoupk.oceansalipg.annual.nc'},
             'T':{'lon':'longitude','var_name':'temp_ym_dpth', 'folder':'time_series', 'name':'xoupk.oceansalipg.annual.nc'}}


# dim = {'S':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'salinity_ym_dpth'},
#        'T':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'temp_ym_dpth'}, 
#        'U':{'lon':'longitude_1', 'lat':'latitude_1', 'depth':'depth_1', 'var':'ucurrTot_ym_dpth'}, 
#        'V':{'lon':'longitude_1', 'lat':'latitude_1', 'depth':'depth_1', 'var':'vcurrTot_ym_dpth'},
#        'Ice':{'lon':'longitude', 'lat':'latitude', 'depth':'surface', 'var':'iceconc_mm_srf'},
#        'Advx':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth'},
#        'Advy':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_1'},
#        'Advz':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_2'}, 
#        'Difx':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_3'},
#        'Dify':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_4'},
#        'Difz':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_5'},
#        'Conv':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_9'},
#        'ML':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_8'},
#        'Fourier':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_11'},
#        'GM':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_14'},
#        'Robert':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_12'}, 
#        'dIce':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_2', 'var':'dSalinitydt_ym_dpth_7'},
#        'Med':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_1', 'var':'dSalinitydt_ym_dpth_13'},
#        'sfc':{'lon':'longitude', 'lat':'latitude', 'depth':'depth_2', 'var':'dSalinitydt_ym_dpth_6'}}

def create_decadal_na(variable):

    print(f"__ Creating decadal time series for {variable}")

    name_input = f"/nfs/see-fs-01_users/eeymr/database/{expt}/{variable['folder']}/{variable['name']}"
    name_output = name_input.replace('annual','decadal').replace(variable['folder'],'decadal').replace('.nc','.na.nc')

    print(f"__ Opening {name_input}")
    ds = xr.open_dataset(name_input)

    print(f"__ Processing")
    ds_out = util.groupby_decade(util.extend_lon(ds[variable['var_name']], n_lon, lon_name=variable['lon']).isel({variable['lon']:slice(na_lon[0], na_lon[1])}))

    print(f"Saving at {name_output}")
    if os.path.exists(name_output):
        os.remove(name_output)
    ds_out.to_netcdf(name_output)


def main():
    # for variable_name in variables.keys():
    for variable_name in ['S']:
        create_decadal_na(variables[variable_name])


main()
