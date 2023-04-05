import numpy as np
import pandas as pd
import xarray as xr 
from divergent_wind import Hadley_Walker

# Load of data and prep
model = 'MOHC/HadGEM3-GC31-LL'
path_to_data = '/badc/cmip6/data/CMIP6/CMIP/'
data_path = path_to_data+model+'/historical/r1i1p1f3/Amon/'
data_dicc={}
for var in ['ua','va']:
	filename = data_path+var+'/gn/latest/'+var+'_Amon_HadGEM3-GC31-LL_historical_r1i1p1f3_gn_195001-201412.nc'
	da = xr.open_dataset(filename)
	data_dicc[var]=da[var]

#Define outputname
outname = 'output/HadGEM3-GC31-LL_historical_'
# Save individual wind components Udiv,Vdiv,Urot,Vrot
save_components=False
# Save  and compute streamfunction
stream_bool=True
# Save and compute omega and mass flux
mass_bool=True

# Actual functions
# Init class as xarray dataarrays U,V
div_wind = Hadley_Walker(data_dicc['ua'],data_dicc['va'],outname)
div_wind.get_components(save_components)
div_wind.compute_diagnostics(stream=stream_bool,mass_omega=mass_bool)
# ERA5 data, use one for each year 

