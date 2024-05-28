import numpy as np
import pandas as pd
from dateutil.relativedelta import relativedelta
import xarray 
from scipy.integrate import cumtrapz
import multiprocessing
import sys,copy
from windspharm.xarray import VectorWind
class Hadley_Walker:
    def __init__(self,dau,dav,outname,multi_flag=False):
        """ Class to compute the divergent wind components and or Walker and Had. Circulation

	    Parameters 
	    ---------- 
	    dau : xarray-DataArray 
		3-dimensional object of xarray for the zonal component of the wind U.  
	    dav : xarray-DataArray 
		3-dimensional object of xarray for the meridional component of the wind V.  
	    outname : string
		Main string to name the output, typically containts the main path and model name.   
	    multi_flag : boolean
		True to indicate user wants to use the multiprocessing routines which may make things faster or also might request too much memory. Default is False.
	 
	    Returns 
	    ------- 
        Class of Hadley_Walker instance 

        """
        # Load file #
        if not isinstance(dau, xarray.DataArray) or not isinstance(dav, xarray.DataArray):
        	raise TypeError('u and v must be xarray.DataArray instances')

        # Define constants
        self.a0 = 6371000
        self.g = 9.81
        self.outname =outname
        self.da_u = dau
        self.da_v = dav

        # Check time coordinate
        self.check_timecoord()

        # Multiprocess data True or False
        self.no_multiprocess=multi_flag

        try:
            self.lambda_lon=dau.longitude
            self.phi_lat=dau.latitude
        except:
            self.lambda_lon=dau.lon
            self.phi_lat=dau.lat
        #self.lat_lon_check()
    def lat_lon_check(self):
        print(self.phi_lat)
        if self.phi_lat.data[0]>self.phi_lat.data[-1]:
            self.phi_lat=np.flip(self.phi_lat)
            array.data=np.flip(array.data,axis=2)
            array.latitude=np.flip(array.latitude,axis=2)
        if self.lambda_lon.data[0]>self.lambda_lon.data[-1]:
            self.lambda_lon=np.flip(self.lambda_lon)
            array.longitude=np.flip(array.longitude)
            array.data=np.flip(array.data,axis=3)
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
            
    def streamfunction(self,array,typo,collapsed=False):
        fact = 1
        if np.max(array.coords[self.pcoord].data)<2000:
            p=array.coords[self.pcoord].data*100
        else:
            p=array.coords[self.pcoord].data
        if p[0]<1:# or p[-1]>=1000:
            print('flipped')
            p=np.flip(p)
            array.data=np.flip(array.data,axis=1)
        else:
            print('not flipped')
        array=array.assign_coords({self.pcoord:p})
        if typo=='hadley':
            fact=2*self.a0*np.pi/9.81
            psi=copy.deepcopy(array)
            psi1=cumtrapz(psi.data,x=p,axis=1,initial=0)
            multip_fact = -1*fact*np.cos(self.phi_lat.data*np.pi/180.)
            psi.data=psi1*multip_fact[:,np.newaxis]

            return psi
        elif typo=='walker':
            fact=-1*self.a0*np.pi/9.81
            psi=copy.deepcopy(array)
            psi1=cumtrapz(psi.data,x=p,axis=1,initial=0)
            psi.data=psi1*fact
            return psi
    def omega_2d(self,coord,array,typo='div'):
        """ Compute Omega and mass flux on each component


        """
        # Correction for pressure if pressure is in hPa 
        # pressure must be in Pa. 
        if np.max(array.coords[self.pcoord].data)<2000:
            print('multiplied p')
            p=array.coords[self.pcoord].data*100
        else:
            p=array.coords[self.pcoord].data
        if p[0]<100 or p[-1]>1e4:
            p=np.flip(p)
            array.data=np.flip(array.data,axis=1)
        else:
            print('not flipped omega')
            #p=np.flip(p)
            #array.data=np.flip(array.data,axis=1)
        print(p)
        array=array.assign_coords({self.pcoord:p})
        cos_alpha = np.cos(np.deg2rad(self.phi_lat)).data
        fact_ag = 1/(self.a0*self.g)
        fact_om = -1/(self.a0*cos_alpha)
        mass_flux = copy.deepcopy(array)
        omega = copy.deepcopy(array)
        if coord=='phi':
            grad_phi = np.gradient(array.data*cos_alpha[:,np.newaxis],np.deg2rad(self.phi_lat).data,axis=2)
            if self.phi_lat.data[0]>self.phi_lat.data[-1]:
                print('reversed phi')
                grad_phi=grad_phi*(-1)
            psi1=cumtrapz(grad_phi,x=p,axis=1,initial=0)
            m_phi = fact_ag*psi1
            mass_flux.data = m_phi
            omega_phi = psi1*fact_om[:,np.newaxis]
            omega.data = omega_phi
            return mass_flux,omega
        elif coord=='lambda':

            grad_lam = np.gradient(array.data,np.deg2rad(self.lambda_lon).data,axis=3)

            psi1=cumtrapz(grad_lam,x=p,axis=1,initial=0)
            m_phi = fact_ag*psi1
            mass_flux.data = m_phi
            omega_phi = psi1*fact_om[:,np.newaxis]
            omega.data = omega_phi

            return mass_flux,omega
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
            xr_psi = self.streamfunction(self.v_div,typo='hadley')
            self.array_add_attrs_Save(xr_psi,'streamfunction',{'long_name':'Meridional streamfunction','units':'kg m**-2 s**-1'},outfil=self.outname+'_stream_phi.nc',save=True)
            psi_lambda = self.streamfunction(self.u_div,typo='walker')
            self.array_add_attrs_Save(psi_lambda,'streamfunction',{'long_name':'Zonal streamfunction','units':'kg m**-2 s**-1'},outfil=self.outname+'_stream_lambda.nc',save=True)

        if mass_omega:
            mass_flux,omega = self.omega_2d('phi',self.v_div)
            omega=self.array_add_attrs_Save(omega,'omega',{'long_name':'Vertical velocity meridional component','units':'Pa s**-1'},outfil=self.outname+'_omega_phi.nc',save=True)
            mass_flux=self.array_add_attrs_Save(mass_flux,'mass_flux',{'long_name':'Meridional mass flux','units':'kg m**-2 s**-1'},outfil=self.outname+'_mass_phi.nc',save=True)
            mass_flux_lambda,omega_lambda = self.omega_2d('lambda',self.u_div)
            omega_lambda=self.array_add_attrs_Save(omega_lambda,'omega',{'long_name':'Vertical velocity zonal component','units':'Pa s**-1'},outfil=self.outname+'_omega_lambda.nc',save=True)
            mass_flux_lambda=self.array_add_attrs_Save(mass_flux_lambda,'mass_flux',{'long_name':'Zonal mass flux','units':'kg m**-2 s**-1'},outfil=self.outname+'_mass_lambda.nc',save=True)

            if returni:
                return mass_flux,omega,mass_flux_lambda,omega_lambda
    def get_uchi(self,it):
        dt = self.da_u.time[it]	
        w = VectorWind(self.da_u[it],self.da_v[it])
        print('decomposing wind time-step ',it)
        udiv, vdiv, urot, vrot = w.helmholtz(truncation=21)
        urot = urot.assign_coords(time=dt)
        vrot = vrot.assign_coords(time=dt)
        urot = urot.expand_dims('time')
        vrot = vrot.expand_dims('time')
        udiv = udiv.assign_coords(time=dt)
        vdiv = vdiv.assign_coords(time=dt)
        udiv = udiv.expand_dims('time')
        vidv = vdiv.expand_dims('time')
        return udiv,vdiv,urot,vrot
    def get_components(self,save):
        M = []
        if self.da_u.isnull().any():
            print('NaN values found, extrapolating')
            self.da_u=self.da_u.interpolate_na(dim='lat',method='linear',fill_value='extrapolate')
            self.da_v=self.da_v.interpolate_na(dim='lat',method='linear',fill_value='extrapolate')
            #if self.da_u.isnull().any():
            #    print(np.where(self.da_u.isnull())[1])
            #    plev = self.da_u.plev
#                self.da_u = self.da_u.where(self.da_u.notnull())
#                self.da_u=self.da_u[plev>100]#interpolate_na(dim='plev',method='linear',fill_value='extrapolate')
#                self.da_v=self.da_v[plev>100]#interpolate_na(dim='plev',method='linear',fill_value='extrapolate')

        if len(self.da_u.shape)>3 and self.da_u.shape[0]==1:
            M=[self.get_uchi(0)]
        elif self.no_multiprocess:
            M=[]
            for it in range(self.da_u.shape[0]):
                M.append(self.get_uchi(0))	
        else:
            with multiprocessing.Pool(6) as pool:
                M = pool.map(self.get_uchi,list(range(self.da_u.shape[0])))

        u_div_list=[]
        v_div_list=[]
        v_rot_list=[]
        u_rot_list=[]
        for result in M:
            udiv,vdiv,urot,vrot=result
            u_div_list.append(udiv)	
            v_div_list.append(vdiv)	
            v_rot_list.append(vrot)	
            u_rot_list.append(urot)	
        full_udv = xarray.concat(u_div_list,dim='time')
        full_vdv = xarray.concat(v_div_list,dim='time')
        full_vrot = xarray.concat(v_rot_list,dim='time')
        full_urot = xarray.concat(u_rot_list,dim='time')
        if save:
            full_udv.to_netcdf(self.outname+'_udiv.nc')
            full_vdv.to_netcdf(self.outname+'_vdiv.nc')
            full_vrot.to_netcdf(self.outname+'_vrot.nc')
            full_urot.to_netcdf(self.outname+'_urot.nc')
        self.u_div=full_udv
        self.v_div=full_vdv
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

