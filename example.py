'''
Example: Using Hadley_Walker to Decompose Circulation

This script demonstrates how to:
  - Load example climate model data (UA, VA components)
  - Prepare climatological monthly means
  - Initialize and configure the `Hadley_Walker` class
  - Compute divergent/rotational components
  - Compute and save Hadley and Walker streamfunctions and mass/Omega diagnostics

Notes:
  - Input DataArrays must have dimensions: (time, pressure, latitude, longitude)
  - Coordinates:
      * time: datetime-like
      * month/clim groupby: using `time.month` for computing monthly climatology
      * latitude: degrees_north (named 'lat' or 'latitude')
      * longitude: degrees_east (named 'lon' or 'longitude')
      * pressure: in hPa or Pa (the class handles unit conversion)
  - Outputs are written to NetCDF files prefixed by `outname`.
'''

import numpy as np
import pandas as pd
import xarray as xr
from decompose_wind import Hadley_Walker

# -----------------------------------------------------------------------------
# 1. Load and inspect example data
# -----------------------------------------------------------------------------
model = 'BCC-ESM1'
path_base = '/badc/cmip6/data/CMIP6/CMIP/'
exp_path = f'{path_base}BCC/{model}/amip/r1i1p1f1/Amon/'

# Variables to load
vars_to_load = ['ua', 'va']  # zonal and meridional wind
clim_dic = {}

for var in vars_to_load:
    # Construct file path
    filename = (
        f"{exp_path}{var}/gn/latest/{var}_Amon_{model}_amip_r1i1p1f1_gn_197901-201412.nc"
    )
    ds = xr.open_dataset(filename)

    # DataArray shape: (time: 36 years x 12 months = 432, level, lat, lon)
    print(f"Loaded {var}: {ds[var].shape}, coords: {list(ds[var].coords)}")

    # Store raw and climatological mean
    clim_dic[var] = ds[var]
    clim_dic[f'{var}_clim'] = ds[var].groupby('time.month').mean(dim='time')
    # Climatology shape: (month=12, level, lat, lon)
    print(f"Climatology for {var}: {clim_dic[f'{var}_clim'].shape}")

# -----------------------------------------------------------------------------
# 2. Configure Hadley_Walker
# -----------------------------------------------------------------------------
# Output file prefix
outname = f'output/climatology_{model}_amip'

# Initialize with climatology arrays
# Note: DataArrays must have dims (month, level, lat, lon) -> treated as 'time'
div_wind = Hadley_Walker(
    dau=clim_dic['ua_clim'],
    dav=clim_dic['va_clim'],
    outname=outname,
    multi_flag=False  # set True to enable parallel processing
)

# -----------------------------------------------------------------------------
# 3. Compute divergent and rotational components
# -----------------------------------------------------------------------------
# This step yields attributes: udiv, vdiv, urot, vrot (all xarray.DataArray)
div_wind.get_components(save=False)
print("udiv shape:", div_wind.udiv.shape)  # dims: (time=12, level, lat, lon)

# -----------------------------------------------------------------------------
# 4. Compute diagnostics: streamfunctions and mass/Omega
# -----------------------------------------------------------------------------
# Flags to control saving and computation
compute_stream = True
compute_mass_omega = True

def compute_diagnostics(dw: Hadley_Walker, stream: bool, mass_omega: bool):
    # Hadley streamfunction (meridional)
    if stream:
        psi_h = dw.streamfunction(dw.vdiv, typo='hadley', collapsed=False)
        psi_w = dw.streamfunction(dw.udiv, typo='walker', collapsed=False)
        psi_h.to_netcdf(f"{dw.outname}_psi_hadley.nc")
        psi_w.to_netcdf(f"{dw.outname}_psi_walker.nc")

    # Mass flux and Omega
    if mass_omega:
        mf_h, om_h = dw.omega_2d('phi', dw.vdiv)
        mf_w, om_w = dw.omega_2d('lambda', dw.udiv)
        mf_h.to_netcdf(f"{dw.outname}_mf_hadley.nc")
        om_h.to_netcdf(f"{dw.outname}_omega_hadley.nc")
        mf_w.to_netcdf(f"{dw.outname}_mf_walker.nc")
        om_w.to_netcdf(f"{dw.outname}_omega_walker.nc")

compute_diagnostics(div_wind, compute_stream, compute_mass_omega)

print("Diagnostics complete. Files written to:")
print(" ", outname + '_psi_hadley.nc', outname + '_psi_walker.nc')
print(" ", outname + '_mf_hadley.nc', outname + '_omega_hadley.nc')
print(" ", outname + '_mf_walker.nc', outname + '_omega_walker.nc')
