import xarray as xr
import numpy as np
import os 

import pylaeoclim_leeds.util_hadcm3 as util
from shapely.geometry import box, Point, Polygon


'''
Calculate the spans for each simulation presented in PhD2
'''

# experiments = ['xoupk', 'tfgbi', 'xouph', 'xoupf', 'xppbf']
experiments = ['xppbf']

print(f"Running {os.path.basename(__file__)}")

# Global variables
database = "/nfs/see-fs-01_users/eeymr/database"
ds_basin = xr.open_dataset("~/work/data/um/basin_hadcm3_glac1d_lgm.nc")

masks_na = {}

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

# Import the AMOC and apply filtering
print("__ Importing and filtering AMOC")
amoc_ts = {}
amocf_ts = {}

for expt in experiments:
    ts = xr.open_dataset(f"/nfs/see-fs-01_users/eeymr/database/{expt}/time_series/{expt}.merid.annual.nc")
    amoc_ts[expt] = ts.Merid_Atlantic.sel(latitude=26.5, method='nearest').max('depth').sortby('t')
    amocf_ts[expt] = filt.process(amoc_ts[expt].values)
    

# Import MLD in key convection sites and apply filtering
print("__ Importing and filtering MLD")
mld = {}
mld_gin, mld_irm = {}, {}
mldf_gin, mldf_irm = {}, {}

for expt in experiments:
    ts = xr.open_dataset(f"{database}/{expt}/perso/{expt}.oceanmixedpf.annual.nc").mixLyrDpth_mm_uo

    if 'unspecified' in ts.dims:
        ts = ts.isel(unspecified=0, drop=True)
    
    mld_gin[expt] = (ts*masks_na['gin']).mean(['longitude','latitude'])
    mld_irm[expt] = (ts*masks_na['irm']).mean(['longitude','latitude'])
    mldf_gin[expt] = filt.process(mld_gin[expt].values)
    mldf_irm[expt] = filt.process(mld_irm[expt].values)


# Define the span boxes in the MLD space
print("__ Calculating spans")

span_regions = {}

span_regions['tfgbi']={'cold':Polygon([(18,30), (27,30), (27,40), (18,40)]), 
                       'merid':Polygon([(70,50), (95,50), (95,90), (75,90), (55,65)]), 
                       'zonal':Polygon([(30,65), (55,65), (75,90), (55,105), (30,105)]),
                       'warming':Polygon([(18,40), (27,40), (70,50), (55,65)]), 
                       'cooling':Polygon([(18,40), (27,40), (55,65), (30,65)])}

span_regions['xppbf']={'cold':Polygon([(18,30), (27,30), (27,40), (18,40)]), 
                       'merid':Polygon([(70,50), (100,50), (100,90), (75,90), (55,65)]), 
                       'zonal':Polygon([(30,65), (55,65), (75,90), (55,105), (30,105)]),
                       'warming':Polygon([(18,40), (27,40), (70,50), (55,65)]), 
                       'cooling':Polygon([(18,40), (27,40), (55,65), (30,65)])}

span_regions['xoupk']={'cold':Polygon([(18,30), (27,30), (27,40), (18,40)]), 
                       'merid':Polygon([(70,50), (100,50), (100,90), (75,90), (55,60)]), 
                       'zonal':Polygon([(30,65), (55,60), (75,90), (55,105), (30,105)]),
                       'warming':Polygon([(18,40), (27,40), (70,50), (55,60)]), 
                       'cooling':Polygon([(18,40), (27,40), (55,60), (30,65)])}

span_regions['xouph']={'cold':box(35,40,45,55), 
                       'warming':box(35,55,52,67), 
                       'merid':box(0,0,0,0), 
                       'zonal':box(35,67,52,95),
                       'cooling':box(35,55,52,67)}

span_regions['xoupf']={'cold':box(10,25,25,55), 
                       'warming':box(25,30,45,55), 
                       'merid':box(45,30,95,55), 
                       'zonal':box(0,0,0,0), 
                       'cooling':box(25,30,45,55)}

spans_conditions = {}


# Define the time steps where each conditions is respected
for expt in experiments:
    
    spans_conditions[expt] = {}
    damocf_ts = np.diff(util.rmean(amocf_ts[expt], 100))
    
    for phase in span_regions[expt].keys():
        
        if phase=='warming':
            spans_conditions[expt][phase] = [ \
                span_regions[expt][phase].contains(Point(mldf_gin[expt][i], mldf_irm[expt][i])) and
                damocf_ts[i]>0
                                       for i in range(len(amocf_ts[expt])-1)]
        elif phase=='cooling':
            spans_conditions[expt][phase] = [ \
                span_regions[expt][phase].contains(Point(mldf_gin[expt][i], mldf_irm[expt][i])) and
                damocf_ts[i]<0
                                       for i in range(len(amocf_ts[expt])-1)]
        else:
            spans_conditions[expt][phase] = [span_regions[expt][phase].contains(
                Point(mldf_gin[expt][i], mldf_irm[expt][i]))
                                       for i in range(len(amocf_ts[expt])-1)]
        
        spans_conditions[expt][phase][:50] = [False]*50 # The first 50 time steps are not allocated
        spans_conditions[expt][phase].append(False) # The last time step is not allocated
        
        
# Because the shorter time series prevent to accurately capture the oscillating filltered time series,
# we crop the result of TFGBI to calculate the phase of XOUPK
# spans_conditions['xoupk'] = {}
# for phase in span_regions['xoupk'].keys(): 
#     spans_conditions['xoupk'][phase] = spans_conditions['tfgbi'][phase][:4000]


# Define the spans (i.e the slice where each condition is respected)
spans = {}

for expt in spans_conditions.keys():
    spans[expt] = {}
    
    for phase in spans_conditions[expt].keys():
        span_stack = []
        for i in range(1,len(spans_conditions[expt][phase])):
            if spans_conditions[expt][phase][i]!=spans_conditions[expt][phase][i-1]:
                span_stack.append(i)
        spans[expt][phase] = [(span_stack[2*i],span_stack[2*i+1]) for i in range(len(span_stack)//2)]


# Saving the spans
print('__ Saving spans')

for expt in experiments:
    sav_folder = f"/nfs/see-fs-01_users/eeymr/database/{expt}/perso/spans"
    util.makedir(sav_folder)

    for phase in span_regions[expt].keys():
        util.np_sav(spans_conditions[expt][phase], f"{sav_folder}/{expt}.span_condition.{phase}")
        util.np_sav(spans[expt][phase], f"{sav_folder}/{expt}.spans.{phase}")
        