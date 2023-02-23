import numpy as np 
import os,glob
import xarray as xr
from divergent_wind import divergent_wind
#path to U and V files
upath = ''
vpath = ''

date='feb10'

year ='2019'
ens ='ens3'

path = '/discover/nobackup/projects/geos_s2s_tc/jorge_work/vars_3d/'+date+'/'+year+'/'+ens+'/'
filelist = np.sort(glob.glob(path+'*nc4'))
outpath='/discover/nobackup/projects/geos_s2s_tc/jorge_work/output/var_files/'+date+'_'+year+'_'+ens
#uda = xr.open_dataarray(upath)
#vda = xr.open_dataarray(vpath)

ds = xr.open_mfdataset(filelist)
uda,vda=ds['U'],ds['V']

# output path
#outpath = ''

# save divergent and rotational wind component output
component_save = False
compute_omega = True
compute_streamfunction = True
# Save omega, mass flux and streamfunction output
save_omega_stream = True
# return diagnostics, False terminates the script
returni=False



#init class
div_wind = divergent_wind(uda[0:3],vda[0:3],outpath)
div_wind.get_components(save=component_save)

ret_list=div_wind.compute_streamfunction(outpath,omega=True,stream=True,returni=False,save=save_omega_stream)

