import numpy as np
import pandas as pd
from dateutil.relativedelta import relativedelta
import xarray 
from scipy.integrate import cumtrapz
import multiprocessing
import sys,copy
from windspharm.xarray import VectorWind
class Hadley_Walker:
    def __init__(self,dau,dav,outname):
        # Load file #
        if not isinstance(dau, xarray.DataArray) or not isinstance(dav, xarray.DataArray):
        	raise TypeError('u and v must be xarray.DataArray instances')
        self.a0 = 6371000
        self.g = 9.81
        self.outname =outname
        self.da_u = dau
        self.pcoord = list(dau.coords)[1]
        print(self.pcoord)
        self.da_v = dav
        array=dau
        try:
            self.lambda_lon=array.longitude
            self.phi_lat=array.latitude
        except:
            self.lambda_lon=array.lon
            self.phi_lat=array.lat

    def streamfunction(self,array,typo,collapsed=False):
        fact = 1
        if np.max(array.coords[self.pcoord].data)<2000:
            p=array.coords[self.pcoord].data*100
        else:
            p=array.coords[self.pcoord].data
            p=np.flip(p)
            array.data=np.flip(array.data,axis=1)
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
        if np.mean(array.coords[self.pcoord].data)<1000:
            p=array.coords[self.pcoord].data*100
        else:
            p=array.coords[self.pcoord].data
            p=np.flip(p)
            array.data=np.flip(array.data,axis=1)

        array=array.assign_coords({self.pcoord:p})

        cos_alpha = np.cos(np.deg2rad(self.phi_lat)).data
        if coord=='phi':
            fact_ag = 1/(self.a0*self.g)
            fact_om = -1/(self.a0*cos_alpha)
            mass_flux = copy.deepcopy(array)
            omega = copy.deepcopy(array)
            if typo=='div':
                grad_phi = np.gradient(array.data*cos_alpha[:,np.newaxis],np.deg2rad(self.phi_lat).data,axis=2)
                psi1=cumtrapz(grad_phi,x=p,axis=1,initial=0)
                m_phi = fact_ag*psi1
                mass_flux.data = m_phi
                omega_phi = psi1*fact_om[:,np.newaxis]
                omega.data = omega_phi

            elif typo=='psi':
                psi = array
                psi1=cumtrapz(psi.data,x=p,axis=1,initial=0)
                psi1=cumtrapz(psi.data,x=p,axis=1,initial=0)

                psi_phi = psi.data*cos_alpha[:,np.newaxis]
                psi_grad = copy.deepcopy(psi)
                psi_grad.data = np.gradient(psi_phi,np.deg2rad(alpha),axis=2)
                psi_grad.data = psi_grad.data*fact[:,np.newaxis]
                psi_mean = psi_grad.collapsed('time',iris.analysis.MEAN)
            return mass_flux,omega
        elif coord=='lambda':
            fact_ag = 1/(self.a0*self.g)
            fact_om = -1/(self.a0*cos_alpha)
            mass_flux = copy.deepcopy(array)
            omega = copy.deepcopy(array)
            if typo=='div':
                grad_phi = np.gradient(array.data,np.deg2rad(self.lambda_lon).data,axis=3)
                psi1=cumtrapz(grad_phi,x=p,axis=1,initial=0)
                m_phi = fact_ag*psi1
                mass_flux.data = m_phi
                omega_phi = psi1*fact_om[:,np.newaxis]
                omega.data = omega_phi

            return mass_flux,omega
    def array_add_attrs_Save(self,da,name,attributes,outfil,save=False):
        da.name=name
        da.attrs=attributes


        if save:
            da.to_netcdf(outfil)
        return da
    def compute_diagnostics(self,stream=False,mass_omega=False,returni=False):
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
    #		uchi, vchi = w.irrotationalcomponent()
        udiv, vdiv, urot, vrot = w.helmholtz(truncation=21)
        #strmf,vp = w.sfvp(truncation=21)
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
        if len(self.da_u.shape)>3 and self.da_u.shape[0]==1:
            M=[self.get_uchi(0)]
        else:
            with multiprocessing.Pool() as pool:
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
        else:
            self.u_div=full_udv
            self.v_div=full_vdv
        return

# Load U, V from ERA5

# Actual functions

# ERA5 data, use one for each year 

