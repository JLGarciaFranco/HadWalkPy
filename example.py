import numpy as np
import pandas as pd
import xarray as xr 
from decompose_wind import Hadley_Walker

# Load of data and prep
model = 'BCC-ESM1'
model_name = 'BCC/BCC-ESM1' 
path_to_data = '/badc/cmip6/data/CMIP6/CMIP/'
data_path = path_to_data+model_name+'/amip/r1i1p1f1/Amon/'
# dictionary to store xarrays
data_dicc={}
for var in ['ua','va']:
	filename = data_path+var+'/gn/latest/'+var+'_Amon_'+model+'_amip_r1i1p1f1_gn_197901-201412.nc'
	da = xr.open_dataset(filename)
	data_dicc[var]=da[var]
	data_dicc[var+'clim']=da[var].groupby('time.month').mean()

# Flags to configure the class that runs the code
# Save individual wind components Udiv,Vdiv,Urot,Vrot
save_components=False
# Save  and compute streamfunction
stream_bool=True
# Save and compute omega and mass flux
mass_bool=True

#Define outputname
outname = 'output/climatology_'+model+'_amip_'

# Actual functions
# create class
div_wind = Hadley_Walker(data_dicc['uaclim'],data_dicc['vaclim'],outname)
# compute wind components
div_wind.get_components(save_components)
# compute diagnostics of the hadley and walker circulations
div_wind.compute_diagnostics(stream=stream_bool,mass_omega=mass_bool)

