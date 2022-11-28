import numpy as np
import matplotlib
matplotlib.use('agg') 
import pandas as pd
from dateutil.relativedelta import relativedelta
import iris
import xarray 
import iris.coord_categorisation
from iris.analysis.calculus import cube_delta
from iris.analysis.calculus import differentiate
import multiprocessing

import sys
sys.path.append('/home/users/jlgarcia/')
from analysis_package import *
from plotting_package import *
from windspharm.xarray import VectorWind
class divergent_wind:
	def __init__(self,dau,dav):
		# Load file #
		self.a0 = 6371000
		self.g = 9.81

		self.da_u = dau
		self.da_v = dav
	def streamfunction(self,array,pcoord,typo,collapsed=False):


		lambda_lon=array.coord('longitude').points
		phi_lat=array.coord('latitude').points
		fact = 1
		if np.mean(array.coord(pcoord).points)<1000:
			p=array.coord(pcoord).points*100
		else:
			p=array.coord(pcoord).points
			p=np.flip(p)
			array.data=np.flip(array.data,axis=1)
		print(p)
		array.coord(pcoord).points=p	

		if typo=='hadley':
			fact=2*self.a0*np.pi/9.81
			psi=copy.deepcopy(array)
			psi1=cumtrapz(psi.data,x=p,axis=1,initial=0)
			multip_fact = -1*fact*np.cos(phi_lat*np.pi/180.)
			psi.data=psi1*multip_fact[:,np.newaxis]

			return psi
		elif typo=='walker':

			fact=-1*self.a0*np.pi/9.81
			psi=copy.deepcopy(array)
			psi1=cumtrapz(psi.data,x=p,axis=1,initial=0)
			psi.data=psi1*fact#multip_fact[:,np.newaxis]
			return psi
			#multip_fact = -1*fact*np.cos(phi_lat*np.pi/180.)
		if collapsed:
			return v.collapsed('time',iris.analysis.MEAN) 
		else:
			return v#.collapsed('time',iris.analysis.MEAN) 
	def omega_2d(self,coord,pcoord,array,typo='div'):
		if np.mean(array.coord(pcoord).points)<1000:
			p=array.coord(pcoord).points*100
		else:
			p=array.coord(pcoord).points
			p=np.flip(p)
			array.data=np.flip(array.data,axis=1)
		array.coord(pcoord).points=p	

		alpha=array.coord('latitude').points
		lambda_lon=array.coord('longitude').points
		cos_alpha = np.cos(np.deg2rad(alpha))
		if coord=='phi':
			fact_ag = 1/(self.a0*self.g)
			fact_om = -1/(self.a0*cos_alpha)
			mass_flux = copy.deepcopy(array)
			omega = copy.deepcopy(array)
			if typo=='div':
				grad_phi = np.gradient(array.data*cos_alpha[:,np.newaxis],np.deg2rad(alpha),axis=2)
				psi1=cumtrapz(grad_phi,x=p,axis=1,initial=0)
				m_phi = fact_ag*psi1
				mass_flux.data = m_phi
				print(np.max(m_phi),np.min(m_phi))	
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
				grad_phi = np.gradient(array.data,np.deg2rad(lambda_lon),axis=3)
				psi1=cumtrapz(grad_phi,x=p,axis=1,initial=0)
				m_phi = fact_ag*psi1
				mass_flux.data = m_phi
				omega_phi = psi1*fact_om[:,np.newaxis]
				omega.data = omega_phi

			return mass_flux,omega
			# Plot check 
			lev = 500*100
			plwindcons=iris.Constraint(pressure_level=lambda cell: cell == lev)
			mass_mean=mass_flux.extract(plwindcons)
			mass_500=mass_mean[5:9].collapsed('time',iris.analysis.MEAN).extract(plwindcons)
			fig = plt.figure(figsize=(12,5))
			cs=cartopy_glbal(fig,1,mass_500,'omega',True,cols=1,rows=1,extent='deeptropics',subtitle=r'$m _{\lambda}$')
			add_map_colorbar('horizontal',[0.1,0.08,0.8,0.009],'kg m s',fig,cs)
			save_fig(fig,'omega_lambda_500','.png')
	def compute_streamfunction(self,outname):
		da = iris.load_cube(outname+'_vdiv.nc')
		print(da)
		psi = self.streamfunction(da,'pressure_level',typo='hadley')
		mass_flux,omega = self.omega_2d('phi','pressure_level',da)

		xr_psi=xr.DataArray.from_iris(psi)
		xr_ome=xr.DataArray.from_iris(omega)
		xr_mas =xr.DataArray.from_iris(mass_flux)

		array_add_attrs_Save(xr_ome,'omega',{'long_name':'Vertical velocity meridional component','units':'Pa s**-1'},outfil=outname+'_omega_phi.nc')
		array_add_attrs_Save(xr_mas,'mass_flux',{'long_name':'Meridional mass flux','units':'kg m**-2 s**-1'},outfil=outname+'_mass_phi.nc')
		array_add_attrs_Save(xr_psi,'streamfunction',{'long_name':'Meridional streamfunction','units':'kg m**-2 s**-1'},outfil=outname+'_stream_phi.nc')

		da = iris.load_cube(outname+'_udiv.nc')
		psi = self.streamfunction(da,'pressure_level',typo='walker')
		mass_flux,omega = self.omega_2d('lambda','pressure_level',da)

		xr_psi_lambda=xr.DataArray.from_iris(psi)
		xr_ome_lambda=xr.DataArray.from_iris(omega)
		xr_mas_lambda =xr.DataArray.from_iris(mass_flux)
		array_add_attrs_Save(xr_ome_lambda,'omega',{'long_name':'Vertical velocity zonal component','units':'Pa s**-1'},outfil=outname+'_omega_lambda.nc')
		array_add_attrs_Save(xr_mas_lambda,'mass_flux',{'long_name':'Zonal mass flux','units':'kg m**-2 s**-1'},outfil=outname+'_mass_lambda.nc')
		array_add_attrs_Save(xr_psi_lambda,'streamfunction',{'long_name':'Zonal streamfunction','units':'kg m**-2 s**-1'},outfil=outname+'_stream_lambda.nc')
	def get_uchi(self,it):
		dt = self.da_u.time[it]	
		w = VectorWind(self.da_u[it],self.da_v[it])
#		uchi, vchi = w.irrotationalcomponent()
		udiv, vdiv, urot, vrot = w.helmholtz(truncation=21)
		#strmf,vp = w.sfvp(truncation=21)
		print(urot)
		urot = urot.assign_coords(time=dt)
		vrot = vrot.assign_coords(time=dt)
		urot = urot.expand_dims('time')
		vrot = vrot.expand_dims('time')
		udiv = udiv.assign_coords(time=dt)
		vdiv = vdiv.assign_coords(time=dt)
		udiv = udiv.expand_dims('time')
		vidv = vdiv.expand_dims('time')
		return udiv,vdiv,urot,vrot
	def get_components(self,outname):
		processes = []
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
		print(M)	
		full_udv = xarray.concat(u_div_list,dim='time')
		full_vdv = xarray.concat(v_div_list,dim='time')
		full_vrot = xarray.concat(v_rot_list,dim='time')
		full_urot = xarray.concat(u_rot_list,dim='time')
		full_udv.to_netcdf(outname+'_udiv.nc')
		full_vdv.to_netcdf(outname+'_vdiv.nc')
		full_vrot.to_netcdf(outname+'_vrot.nc')
		full_urot.to_netcdf(outname+'_urot.nc')
		return
		la=newvar.coords['latitude']
		arry=newvar.loc[dict(latitude=la[(la > -32) & (la < 32)])]
		marry = arry.mean('time').sel(level=500)
		print(marry)
		fig=plt.figure(figsize=(12,8))
		marry.plot.contourf(cmap='RdBu_r')
		save_fig(fig,'div_500','.png')
#		eq_avg = arry.mean(['latitude'])

# Load U, V from ERA5
dataset_name = 'era5'
u_cube=load_ecmwf('era5monthly','x_wind')
v_cube=load_ecmwf('era5monthly','y_wind')
year = int(sys.argv[1])
outname = '/gws/nopw/j04/aopp/jorgelgf/aux_cubes/div_wind/'+dataset_name+'_'+str(year)
newt_u = gtvec(u_cube)
newt_v = gtvec(v_cube)

da_u = xarray.DataArray.from_iris(u_cube[newt_u.year==year])
da_v = xarray.DataArray.from_iris(v_cube[newt_v.year==year])


# Actual functions
div_wind = divergent_wind(da_u,da_v)
div_wind.get_components(outname)
div_wind.compute_streamfunction(outname)
# ERA5 data, use one for each year 

