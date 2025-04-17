'''
Module: circulation_split.py

Decompose atmospheric circulation into Hadley and Walker components
using the Schwendike et al. (2014) methodology.

This module provides a `CirculationDecomposer` class that:
  - Loads zonal and meridional wind fields
  - Computes the zonal-mean streamfunction
  - Performs decomposition to isolate Hadley and Walker circulations
  - Supports multiprocessing for batch processing of multiple time snapshots

References
----------
Schwendike, J., Shepherd, T. G., & Polichtchouk, I. (2014). Separation of the Walker circulation from the Hadley circulation
'''
import numpy as np
import pandas as pd
from dateutil.relativedelta import relativedelta
import xarray as xr
from scipy.integrate import cumtrapz
import multiprocessing
import sys,copy
from windspharm.xarray import VectorWind
class Hadley_Walker:
    """
    Class to compute divergent wind components and split the Walker and Hadley circulations.

    Parameters
    ----------
    dau : xr.DataArray
        3D DataArray of the zonal wind component (dimensions: time, lat, lon).
    dav : xr.DataArray
        3D DataArray of the meridional wind component (dimensions: time, lat, lon).
    outname : str
        Base path and model identifier for output files.
    multi_flag : bool, optional
        If True, enable multiprocessing routines for faster computation (may increase memory usage).
        Defaults to False.

    Attributes
    ----------
    a0 : float
        Earth radius in meters (6371000 m).
    g : float
        Gravitational acceleration (9.81 m/s²).
    outname : str
        Output filename base.
    da_u : xr.DataArray
        Zonal wind input.
    da_v : xr.DataArray
        Meridional wind input.
    no_multiprocess : bool
        Flag indicating whether multiprocessing is disabled.
    lambda_lon : xr.DataArray
        Longitude coordinate.
    phi_lat : xr.DataArray
        Latitude coordinate.

    Methods
    -------
    check_timecoord()
        Verify and standardize the time coordinate in input DataArrays.
    compute_divergent_components()
        Calculate divergent part of the wind field.
    split_circulations()
        Perform the Hadley–Walker decomposition.
    run_parallel(times)
        Execute decomposition over multiple time slices in parallel.
    save_results()
        Save streamfunctions or metrics to NetCDF or text.
    """

    def __init__(
        self,
        dau: xr.DataArray,
        dav: xr.DataArray,
        outname: str,
        multi_flag: bool = False
    ):
        # Validate inputs
        if not isinstance(dau, xr.DataArray) or not isinstance(dav, xr.DataArray):
            raise TypeError("`dau` and `dav` must be xarray.DataArray instances")

        # Fundamental constants
        self.a0 = 6371000.0
        self.g = 9.81
        self.outname = outname

        # Assign wind fields
        self.da_u = dau
        self.da_v = dav

        # Check and standardize time coordinate
        self.check_timecoord()

        # Multiprocessing flag (False disables parallel routines)
        self.no_multiprocess = not multi_flag

        # Extract spatial coordinates
        try:
            self.lambda_lon = dau.longitude
            self.phi_lat = dau.latitude
        except AttributeError:
            # fallback to alternative names
            self.lambda_lon = dau.lon
            self.phi_lat = dau.lat

    def check_timecoord(self):
        """Check the time coordinate in case the input are monthly or day of year climatological means
        
        """
        lat, lat_dim = self._find_latitude_coordinate(self.da_u)
        lon, lon_dim = self._find_longitude_coordinate(self.da_u)
        plev_c, plev_dim = self._find_pressure_coordinate(self.da_u)
        self.pcoord = plev_c.name#list(self.da_u.coords)[1] 
        print(self.pcoord)
        print(lat_dim)
        if 'time' in self.da_u.coords.keys() and 'time' in self.da_v.coords.keys():
            self.month_flag = False
            return 
        else:
            if 'month' in self.da_u.coords.keys() and 'month' in self.da_v.coords.keys():
            # Rename month to time 
                self.da_u=self.da_u.rename({'month':'time'})
                self.da_v=self.da_v.rename({'month':'time'})
            elif 'dayofyear' in self.da_u.coords.keys(): 
                self.da_u=self.da_u.rename({'dayofyear':'time'})
                self.da_v=self.da_v.rename({'dayofyear':'time'})

            # Find the current order 

            order = list(range(self.da_u.ndim))
            order.remove(lat_dim)
            order.remove(plev_dim)
            order.remove(lon_dim)
            order.insert(1, plev_dim)
            order.insert(2, lat_dim)
            order.insert(3, lon_dim)
            # Reorder 
            reorder = [order.index(i) for i in range(self.da_u.ndim)]
            apiorder = [self.da_u.dims[i] for i in order]
            self.da_u = self.da_u.transpose(*apiorder)	
            self.da_v = self.da_v.transpose(*apiorder)	

            # Reindex
            current_indexes = self.da_u.indexes
            reordered_indexes = {index_name: current_indexes[index_name] for index_name in apiorder}
            self.da_u = self.da_u.reindex(reordered_indexes)
            self.da_v = self.da_v.reindex(reordered_indexes)  
            self.month_flag = True

            return
    def streamfunction(
        self,
        array: xr.DataArray,
        typo: str,
        collapsed: bool = False
    ) -> xr.DataArray:
        """
        Compute mass streamfunction for Hadley or Walker component by vertical integration.
        """
        p = array.coords[self.pcoord].data
        if np.max(p)<2000:
            p = p*100.0
        if p[0]<1:
            p = np.flip(p)
            array = array.sortby(self.pcoord, ascending=False)
        array = array.assign_coords({self.pcoord: p})
        # Vertical integration
        psi1 = cumtrapz(array.data, x=p, axis=1, initial=0)
        if typo=='hadley':
            fact = 2*self.a0*np.pi/self.g
            cos_lat = np.cos(np.deg2rad(self.phi_lat.data))
            # Broadcast to (time, plev, lat, lon)
            cos4d = cos_lat[np.newaxis, np.newaxis, :, np.newaxis]
            psi_array = psi1 * (-fact) * cos4d
            psi = xr.DataArray(psi_array, coords=array.coords, dims=array.dims)
        elif typo=='walker':
            fact = -self.a0*np.pi/self.g
            psi_array = psi1 * fact
            psi = xr.DataArray(psi_array, coords=array.coords, dims=array.dims)
        else:
            raise ValueError("typo must be 'hadley' or 'walker'")
        if collapsed:
            return psi.mean(dim=('lat','lon'))
        return psi

    def omega_2d(
        self,
        coord: str,
        array: xr.DataArray,
        typo: str = 'div'
    ) -> tuple:
        """
        Compute vertical mass flux and Omega for each circulation component.

        Parameters
        ----------
        coord : {'phi', 'lambda'}
            Direction of gradient: 'phi' for meridional, 'lambda' for zonal.
        array : xr.DataArray
            Divergent or rotational wind component (dimensions: time, plev, lat, lon).
        typo : str, optional
            Component label, unused in calculation. Default is 'div'.

        Returns
        -------
        mass_flux : xr.DataArray
            Vertical mass flux for the component.
        omega : xr.DataArray
            Omega (vertical velocity) for the component.
        """
        # Ensure pressure coordinate in Pa and correct ordering
        p = array.coords[self.pcoord].data
        if np.max(p) < 2000:
            p = p * 100.0
        if p[0] < 100 or p[-1] > 1e4:
            p = np.flip(p)
            array = array.sortby(self.pcoord, ascending=False)
        array = array.assign_coords({self.pcoord: p})

        # Prepare latitude and longitude arrays
        phi = self.phi_lat.values
        lam = self.lambda_lon.values
        cos_phi = np.cos(np.deg2rad(phi))
        # Broadcast for full 4D multiplication
        cos4d = cos_phi[np.newaxis, np.newaxis, :, np.newaxis]

        # Coefficients
        fact_ag = 1.0 / (self.a0 * self.g)
        # Broadcast cos4d in denominator for Omega
        fact_om = -1.0 / (self.a0 * cos4d)

        # Initialize outputs
        mass_flux = array.copy(deep=True)
        omega = array.copy(deep=True)

        # Compute gradient
        if coord == 'phi':
            data_mod = array.data * cos4d
            grad = np.gradient(data_mod, np.deg2rad(phi), axis=2)
        elif coord == 'lambda':
            grad = np.gradient(array.data, np.deg2rad(lam), axis=3)
        else:
            raise ValueError("coord must be 'phi' or 'lambda'")

        # Vertical integration
        psi1 = cumtrapz(grad, x=p, axis=1, initial=0)

        # Build outputs
        mass_flux.data = fact_ag * psi1
        omega.data = psi1 * fact_om

        return mass_flux, omega 
    def array_add_attrs_Save(self,da,name,attributes,outfil,save=False):
        """ Saving output with specified attributes

	    Parameters 
	    ---------- 
	    da : xarray-DataArray 
		xarray.DataArray to be saved  
	    name : string 
		Name of field.  
	    attributes : dict
		Dictionary of xarray attributes.   

	    multi_flag : boolean
		True to indicate user wants to use the multiprocessing routines which may make things faster or also might request too much memory. Default is False.
	 
	    Returns 
	    ------- 

        """
        da.name=name
        da.attrs=attributes
#        da[self.pcoord]=da[self.pcoord]/100.
        if save:
            da.to_netcdf(outfil)
        return da
    def compute_diagnostics(self,stream=False,mass_omega=False,returni=False):
        print('computing streamfunction and omega ')
        if stream:
            xr_psi = self.streamfunction(self.vdiv,typo='hadley')
            self.array_add_attrs_Save(xr_psi,'streamfunction',{'long_name':'Meridional streamfunction','units':'kg m**-2 s**-1'},outfil=self.outname+'_stream_phi.nc',save=True)
            psi_lambda = self.streamfunction(self.udiv,typo='walker')
            self.array_add_attrs_Save(psi_lambda,'streamfunction',{'long_name':'Zonal streamfunction','units':'kg m**-2 s**-1'},outfil=self.outname+'_stream_lambda.nc',save=True)

        if mass_omega:
            mass_flux,omega = self.omega_2d('phi',self.vdiv)
            omega=self.array_add_attrs_Save(omega,'omega',{'long_name':'Vertical velocity meridional component','units':'Pa s**-1'},outfil=self.outname+'_omega_phi.nc',save=True)
            mass_flux=self.array_add_attrs_Save(mass_flux,'mass_flux',{'long_name':'Meridional mass flux','units':'kg m**-2 s**-1'},outfil=self.outname+'_mass_phi.nc',save=True)
            mass_flux_lambda,omega_lambda = self.omega_2d('lambda',self.udiv)
            omega_lambda=self.array_add_attrs_Save(omega_lambda,'omega',{'long_name':'Vertical velocity zonal component','units':'Pa s**-1'},outfil=self.outname+'_omega_lambda.nc',save=True)
            mass_flux_lambda=self.array_add_attrs_Save(mass_flux_lambda,'mass_flux',{'long_name':'Zonal mass flux','units':'kg m**-2 s**-1'},outfil=self.outname+'_mass_lambda.nc',save=True)

            if returni:
                return mass_flux,omega,mass_flux_lambda,omega_lambda
    def get_uchi(self, it: int):
        """
        Perform Helmholtz decomposition of the wind at a given time index.

        Parameters
        ----------
        it : int
            Index along the time dimension for which to compute decomposition.

        Returns
        -------
        udiv : xr.DataArray
            Divergent component of zonal wind at time index `it`.
        vdiv : xr.DataArray
            Divergent component of meridional wind at time index `it`.
        urot : xr.DataArray
            Rotational (non-divergent) component of zonal wind at time index `it`.
        vrot : xr.DataArray
            Rotational (non-divergent) component of meridional wind at time index `it`.

        Notes
        -----
        Uses `windspharm.xarray.VectorWind.helmholtz` with default spectral truncation of 21.
        The returned fields are assigned the original time coordinate and expanded to 3D.
        """
        # Extract timestamp
        dt = self.da_u.time[it]

        # Initialize VectorWind for this timestep
        w = VectorWind(self.da_u.isel(time=it), self.da_v.isel(time=it))
        print(f"Decomposing wind at time index {it}, timestamp {dt.values}")

        # Helmholtz decomposition: (udiv, vdiv, urot, vrot)
        udiv, vdiv, urot, vrot = w.helmholtz(truncation=21)

        # Attach time coordinate and restore time dimension
        def _assign_time(arr):
            arr = arr.assign_coords(time=dt)
            return arr.expand_dims('time')

        udiv = _assign_time(udiv)
        vdiv = _assign_time(vdiv)
        urot = _assign_time(urot)
        vrot = _assign_time(vrot)

        return udiv, vdiv, urot, vrot
    def get_components(self, save: bool = False) -> None:
        """
        Compute divergent and rotational components for all time steps,
        optionally saving each component to NetCDF.

        Parameters
        ----------
        save : bool, optional
            If True, write full udiv, vdiv, urot, and vrot arrays to NetCDF
            files named `<outname>_udiv.nc`, etc. Default is False.

        Notes
        -----
        - Handles NaNs by latitudinal interpolation if present.
        - Uses multiprocessing if `multi_flag` was set in constructor.
        """
        # Handle NaNs by interpolating along latitude
        if self.da_u.isnull().any():
            print('NaN values found, performing latitudinal interpolation')
            self.da_u = self.da_u.interpolate_na(dim='lat', method='linear', fill_value='extrapolate')
            self.da_v = self.da_v.interpolate_na(dim='lat', method='linear', fill_value='extrapolate')

        # List to collect results
        M = []
        nt = self.da_u.sizes['time']

        # Single timestamp case
        if nt == 1:
            M = [self.get_uchi(0)]
        # Serial loop if multiprocessing disabled
        elif self.no_multiprocess:
            for it in range(nt):
                M.append(self.get_uchi(it))
        # Parallel execution
        else:
            with mp.Pool() as pool:
                M = pool.map(self.get_uchi, list(range(nt)))

        # Separate and concatenate results
        u_div_list, v_div_list, u_rot_list, v_rot_list = zip(*M)
        full_udv = xr.concat(u_div_list, dim='time')
        full_vdv = xr.concat(v_div_list, dim='time')
        full_urot = xr.concat(u_rot_list, dim='time')
        full_vrot = xr.concat(v_rot_list, dim='time')

        # Save to NetCDF if requested
        if save:
            full_udv.to_netcdf(f"{self.outname}_udiv.nc")
            full_vdv.to_netcdf(f"{self.outname}_vdiv.nc")
            full_urot.to_netcdf(f"{self.outname}_urot.nc")
            full_vrot.to_netcdf(f"{self.outname}_vrot.nc")

        # Store results on instance
        self.udiv = full_udv
        self.vdiv = full_vdv
        self.urot = full_urot
        self.vrot = full_vrot
        return
	
    def _find_latitude_coordinate(self,array):
        return self._find_coord_and_dim(
        	array,
	        lambda c: (c.name in ('latitude', 'lat') or
		   c.attrs.get('units') == 'degrees_north' or
		   c.attrs.get('axis') == 'Y'),
	         'latitude')


    def _find_coord_and_dim(self,array, predicate, name):
        """
        Find a dimension coordinate in an `xarray.DataArray` that satisfies
        a predicate function.

        """
        candidates = [coord
    		for coord in [array.coords[n] for n in array.dims]
    			if predicate(coord)]

        if not candidates:
    	    raise ValueError('cannot find a {!s} coordinate'.format(name))
        if len(candidates) > 1:
            msg = 'multiple {!s} coordinates are not allowed'
            raise ValueError(msg.format(name))
        coord = candidates[0]
        dim = array.dims.index(coord.name)
        return coord, dim
    def _find_longitude_coordinate(self,array):
        """Find a longitude dimension coordinate in an `xarray.DataArray`."""
        return self._find_coord_and_dim(
                array,
                lambda c: (c.name in ('longitude', 'lon') or
                c.attrs.get('units') == 'degrees_east' or
                c.attrs.get('axis') == 'X'),
                'longitude')
    def _find_pressure_coordinate(self,array):
        """Find a longitude dimension coordinate in an `xarray.DataArray`."""
        return self._find_coord_and_dim(
	    array,
		lambda c: (c.name in ('plev', 'level','lev','pressure','air_pressure') or
		   c.attrs.get('units') == 'hPa' or
		   c.attrs.get('axis') == 'Z'),
	    'pressure')

