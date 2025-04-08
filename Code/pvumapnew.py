import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import cartopy.crs as ccrs
import cartopy
from datetime import datetime,timedelta
import copy
import math
import sys
import xarray as xr
from geopy import distance
from matplotlib.path import Path
from scipy import ndimage 
from shapely.geometry import Point, Polygon
import glob
import cftime
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset 
from scipy.stats import binned_statistic
from scipy.interpolate import interp1d
import nc_time_axis
import os
os.environ['ESMFMKFILE']='/glade/work/ajacobs/conda-envs/asp2023/lib/esmf.mk'
import xesmf as xe
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
import weather_modules
from scipy.interpolate import griddata
from RegriddingFunctions import hs2p

fv1=xr.open_dataset('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h1.2020-07-28-00000.nc')
fv2=xr.open_dataset('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h2.2020-07-28-00000.nc')
fv3=xr.open_dataset('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h3.2020-07-28-00000.nc')
fv_lin_centers=xr.open_dataset('all_nudge_abl_cyclone_centers.nc')
'''
2016-08-27 00:00:00 cyclone_1764
2018-08-24 00:00:00 cyclone_2765
2020-05-12 18:00:00 cyclone_3566
2013-09-15 18:00:00 cyclone_361
2021-09-22 06:00:00 cyclone_4220
2022-09-17 00:00:00 cyclone_4709
2022-09-20 18:00:00 cyclone_4716
2022-09-25 06:00:00 cyclone_4724
2014-09-26 18:00:00 cyclone_843
2014-09-29 18:00:00 cyclone_845

2013-07-24 06:00:00 cyclone_288
2014-05-16 00:00:00 cyclone_642
2016-06-08 12:00:00 cyclone_1666
2020-05-09 06:00:00 cyclone_3562
2020-05-12 18:00:00 cyclone_3566
2020-07-28 06:00:00 cyclone_3672
'''
#fv1=xr.open_dataset('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_May_Sept2022_fv09/linear_full_nudge/atm/hist/fv9.lin.cam.h1.nc').load()
#fv2=xr.open_dataset('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_May_Sept2022_fv09/linear_full_nudge/atm/hist/fv9.lin.cam.h2.nc').load()
#fv3=xr.open_dataset('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_May_Sept2022_fv09/linear_full_nudge/atm/hist/fv9.lin.cam.h3.nc').load()
def p_interp(pres,var):
#and np.nanmax(theta)>=850:
    f=interp1d(pres,var,bounds_error=False,
                kind='linear',
                fill_value=(np.nan))
    x=fv1.ilev
    return f(x)
def pv_interp(pv,z):
    f=interp1d(pv,z,bounds_error=False,
                kind='linear',
                fill_value=(np.nan))
    return f(2)
def findDistance(lat1,lon1,lat2,lon2):
        center = (lat1,lon1)
        point = (lat2,lon2)
        dist = distance.great_circle(center,point).km
        return dist
pressure=(fv1.hyam*fv1.P0+fv1.hybm*fv2.PS)/100
T= fv1.T
U=fv1.U
V=fv1.V
Q=fv1.Q
Z3=fv3.Z3
ilev=fv1.ilev
lat=fv1.lat.values[130:]
lon=fv1.lon.values
PSL=fv2.PSL
time=fv2.time
fv1.close()
fv2.close()
fv3.close()
for a in range (0,24,6):
    campressl=pressure.transpose("time", "lev", "lat", "lon") # hybrid coord press levels
    theta=weather_modules.temp_to_theta(T, campressl*100)
    theta_pres=hs2p(theta[a:a+1,:,130:,:].data,campressl[a:a+1,:,130:,:].data,ilev.data)
    u_pres=hs2p(U[a:a+1,:,130:,:].data,campressl[a:a+1,:,130:,:].data,ilev.data)
    v_pres=hs2p(V[a:a+1,:,130:,:].data,campressl[a:a+1,:,130:,:].data,ilev.data)
    t_pres=hs2p(T[a:a+1,:,130:,:].data,campressl[a:a+1,:,130:,:].data,ilev.data)
    q_pres=hs2p(Q[a:a+1,:,130:,:].data,campressl[a:a+1,:,130:,:].data,ilev.data)
    z_pres=hs2p(Z3[a:a+1,:,130:,:].data,campressl[a:a+1,:,130:,:].data,ilev.data)
    plev=np.array([300])
    z_pres_300=hs2p(fv3.Z3[a:a+1,:,130:,:].data,campressl[a:a+1,:,130:,:].data,plev)/10
    press_levels3d = np.tile(ilev.data[:, np.newaxis, np.newaxis], (1, 62, 288))
    epv=weather_modules.epv_sphere(theta_pres,press_levels3d*100,u_pres,v_pres,lat,lon)
    pvu_theta=weather_modules.interp2pv(epv, theta_pres, 2)
    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot( projection=ccrs.NorthPolarStereo())
    cs=ax1.pcolormesh(lon,lat[:-2],pvu_theta[:-2],cmap='plasma',vmin=260,vmax=360,transform=ccrs.PlateCarree())
    #cs=ax1.pcolormesh(lon,lat[:-2],z_pres_300[:-2,:],cmap='plasma',vmin=850,vmax=1000,transform=ccrs.PlateCarree())
    cs1=ax1.contour(lon,lat[:-2],PSL[a,130:-2,:]/100,colors='k',levels=np.arange(950,1020,8),transform=ccrs.PlateCarree())
    plt.clabel(cs1, inline_spacing=0, fontsize=10, fmt="%.0f")
    ax1.scatter(fv_lin_centers.lon.values[np.where(np.logical_and(fv_lin_centers.datetime.values==time[a].values,fv_lin_centers.label.values=='cyclone_3672'))[0]],fv_lin_centers.lat.values[np.where(np.logical_and(fv_lin_centers.datetime.values==time[a].values,fv_lin_centers.label.values=='cyclone_3672'))[0]],transform=ccrs.PlateCarree(),color='cyan',marker='x')
    plt.colorbar(cs,fraction=0.046, pad=0.04,label='[K]',extend='both')
    ax1.set_extent([-180,180,55,90], crs=ccrs.PlateCarree()) 
    ax1.coastlines()
    ax1.gridlines()
    title='2PVU Potential Temperature and Sea Level Pressure'+'\n'+'time='+str(time[a].values)
    #title='Geopotential Height at 300 hPa and Sea Level Pressure'+'\n'+'time='+str(time[a].values)
    plt.title(title)
    save_fn =  '%s_[%s].png' % ("2pvu_cyc006_25",str(int((a/6))).zfill(3))
    plt.savefig(save_fn)
