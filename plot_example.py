'''
plot_mass_flux_simple.py

Simplified visualization of Hadley/Walker mass flux without custom dependencies.

This script:
  - Loads precomputed mass flux NetCDF files (`mass_phi`, `mass_lambda`).
  - Filters by DJF (Dec-Jan-Feb) and JJA (Jun-Jul-Aug).
  - Computes time-mean at 50000 Pa.
  - Plots a 2×2 map grid with Cartopy and Matplotlib only.

Requirements:
  - xarray
  - numpy
  - matplotlib
  - cartopy

Usage:
  python plot_mass_flux_simple.py
'''

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# -------------------- Configuration --------------------
model = 'BCC-ESM1'
experiment = 'amip'
outprefix = f'output/climatology_{model}_{experiment}'

# Variables and labels
variables = ['mass_phi', 'mass_lambda']
var_labels = ['Meridional mass flux ($M_\phi$)',
              'Zonal mass flux ($M_\lambda$)']

# Seasons and months
seasons = ['DJF', 'JJA']
season_months = {'DJF': [12, 1, 2], 'JJA': [6, 7, 8]}

# Plot settings
plev_plot = 50000  # Pa
levels = np.arange(-7, 7.1, 0.5)
cmap = plt.get_cmap('RdBu_r')
scale_factors = {'mass_phi': -1000, 'mass_lambda': 1000}

# -------------------- Plotting --------------------
fig = plt.figure(figsize=(12, 6), dpi=150)

for i, season in enumerate(seasons):
    months = season_months[season]
    for j, var in enumerate(variables):
        # subplot index
        ax = fig.add_subplot(2, 2, i*2 + j + 1,
                             projection=ccrs.PlateCarree())

        # Load DataArray
        fname = f'{outprefix}__{var}.nc'
        da = xr.open_dataarray(fname)
        print(da)
        # Ensure datetime time coordinate
        if not np.issubdtype(da.time.dtype, np.datetime64):
            # Select seasonal months
            da_season = da.sel(time=da.time.isin(months))
            #raise ValueError('`time` coordinate must be datetime')
        else:
            da_season = da.sel(time=da.time.dt.month.isin(months))



        # Mean at pressure level
        da_mean = da_season.sel(plev=plev_plot, method='nearest')
        da_mean = da_mean.mean(dim='time', skipna=True)

        # Apply scale
        da_plot = da_mean * scale_factors[var]

        # Plot
        mesh = ax.contourf(
            da_plot.lon, da_plot.lat, da_plot.values,
            levels=levels, cmap=cmap, extend='both',
            transform=ccrs.PlateCarree()
        )
        ax.coastlines(linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3)
        ax.set_title(f'{model} {season}\n{var_labels[j]}', fontsize=10)

# Colorbar
cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.02])
cbar = fig.colorbar(mesh, cax=cbar_ax, orientation='horizontal')
cbar.set_label('[10$^{-3}$ kg m$^{-2}$ s$^{-1}$]')

plt.tight_layout(rect=[0, 0.07, 1, 1])

# Save figure
outfile = f'plots/{model}_{experiment}_mass_flux_simple.png'
plt.savefig(outfile, bbox_inches='tight')
print(f'Figure saved to {outfile}')
