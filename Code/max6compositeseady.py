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
from scipy.interpolate import griddata,Rbf
import nc_time_axis
import os
os.environ['ESMFMKFILE']='/glade/work/ajacobs/conda-envs/asp2023/lib/esmf.mk'
import xesmf as xe
import pandas as pd
import seaborn as sns
import warnings
import metpy.calc as mpcalc
from metpy.units import units
import weather_modules
from RegriddingFunctions import hs2p
warnings.filterwarnings("ignore")
warm_core=[]
cold_core=[]
def findDistance(lat1,lon1,lat2,lon2):
        center = (lat1,lon1)
        point = (lat2,lon2)
        dist = distance.great_circle(center,point).km
        return dist
fv_lin_centers_2=xr.open_dataset('all_nudge_abl_cyclone_centers.nc')
endd=int(fv_lin_centers_2.label.values[-1][-4:])
fv_lin_centers_2=fv_lin_centers_2.sortby(['label'])
for a in [3562]: #3566,3562,288,642,1666,3672
    label = 'cyclone_'+str(a).zfill(3)
    times= fv_lin_centers_2.datetime.values[np.where(fv_lin_centers_2.label.values==label)[0]]
    eady_500_fa=[]
    fv_ind=[]
    for time in times:
        year=str(time.year)
        month=str(time.month).zfill(2)
        day=str(time.day).zfill(2)
        file='/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h2.{y}-{m}-{d}-00000.nc'.format(y=year,m=month,d=day)
        fv2=xr.open_dataset(file)
        file1='/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h1.{y}-{m}-{d}-00000.nc'.format(y=year,m=month,d=day)
        fv1=xr.open_dataset(file1)
        file2='/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h3.{y}-{m}-{d}-00000.nc'.format(y=year,m=month,d=day)
        fv3=xr.open_dataset(file2)
        b=int(time.hour)
        pressure=(fv1.hyam*fv1.P0+fv1.hybm*fv2.PS[b:b+1])/100
        T= fv1.T[b:b+1]
        U=fv1.U[b:b+1]
        V=fv1.V[b:b+1]
        Q=fv1.Q[b:b+1]
        Z3=fv3.Z3[b:b+1]
        ilev=fv1.ilev
        lat=fv1.lat.values[130:]
        lon=fv1.lon.values
        PSL=fv2.PSL[b,130:,:]

        plev=np.array([1000,925,850,700,600,500,400,300,250,100])
        campressl=pressure.transpose("time", "lev", "lat", "lon") # hybrid coord press levels
        theta=weather_modules.temp_to_theta(T, campressl*100)
        theta_pres=hs2p(theta[:,:,130:,:].data,campressl[:,:,130:,:].data,plev)
        u_pres=hs2p(U[:,:,130:,:].data,campressl[:,:,130:,:].data,plev)
        z_pres=hs2p(Z3[:,:,130:,:].data,campressl[:,:,130:,:].data,plev)
        del(campressl)
        clat= fv_lin_centers_2.lat.values[np.where(np.logical_and(fv_lin_centers_2.label.values==label,fv_lin_centers_2.datetime.values==time))[0][0]]
        clon= fv_lin_centers_2.lon.values[np.where(np.logical_and(fv_lin_centers_2.label.values==label,fv_lin_centers_2.datetime.values==time))[0][0]]
        slp_section=[]
        indices=[]
        for count, s in enumerate(PSL):
            for c2, t in enumerate(s):
                dist = findDistance(clat,clon,t.lat,t.lon)
     #   print(dist)
                if dist <=500:
                    slp_section.append(t)
                    indices.append([count,c2])
        indices=np.array(indices)    
        N=np.sqrt(-9.81/theta_pres[:,indices[:,0],indices[:,1]].squeeze()*np.gradient(theta_pres[:,indices[:,0],indices[:,1]].squeeze(),plev,axis=0))
        f=mpcalc.coriolis_parameter(np.deg2rad(lat[indices[:,0]])).magnitude
        f,Y=np.meshgrid(f,np.zeros(len(indices[:,1])))
        dudz=np.gradient(u_pres[:,indices[:,0],indices[:,1]].squeeze(),plev,axis=0)
        eady_500_fa.append(np.mean(0.3098*abs(f)*abs(dudz[5])/N[5]))
        fv_ind.append(time)
   
    
    time=np.array(fv_ind)[np.where(np.array(eady_500_fa)==np.max(np.array(eady_500_fa)))[0][0]]
    year=str(time.year)
    month=str(time.month).zfill(2)
    day=str(time.day).zfill(2)
    if int(month) not in [5,6,7,8,9]:
        continue
    file='/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h2.{y}-{m}-{d}-00000.nc'.format(y=year,m=month,d=day)
    fv2=xr.open_dataset(file)
    file1='/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h1.{y}-{m}-{d}-00000.nc'.format(y=year,m=month,d=day)
    fv1=xr.open_dataset(file1)
    file2='/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h3.{y}-{m}-{d}-00000.nc'.format(y=year,m=month,d=day)
    fv3=xr.open_dataset(file2)
    b=hour=int(time.hour)
    pressure=(fv1.hyam*fv1.P0+fv1.hybm*fv2.PS)/100
    campressl=pressure.transpose("time", "lev", "lat", "lon") # hybrid coord press levels
    ilev=fv1.ilev
    lat=fv1.lat.values[130:]
    lon=fv1.lon.values
    theta=weather_modules.temp_to_theta(fv1.T, campressl*100)
    theta_pres=hs2p(theta[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,ilev.data)
    u_pres=hs2p(fv1.U[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,ilev.data)
    v_pres=hs2p(fv1.V[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,ilev.data)
    t_pres=hs2p(fv1.T[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,ilev.data)
    q_pres=hs2p(fv1.Q[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,ilev.data)
    #z_pres=hs2p(fv3.Z3[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,ilev.data)
    press_levels3d = np.tile(ilev.data[:, np.newaxis, np.newaxis], (1, 62, 288))
    plev=np.array([300])
    #z_pres_300=hs2p(fv3.Z3[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,plev)/10
    epv=weather_modules.epv_sphere(theta_pres,press_levels3d*100,u_pres,v_pres,lat,lon)
    pvu_theta=weather_modules.interp2pv(epv, theta_pres, 2)
    plev=np.array([850])
    u_pres_850=hs2p(fv1.U[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,plev)
    v_pres_850=hs2p(fv1.V[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,plev)
    del(campressl)
    u850=np.sqrt(u_pres_850**2+v_pres_850**2)
    t=fv2.T850[hour,130:,:]
    slp=(fv2.T850)[hour,130:,:]
    tpw=fv2.TMQ[hour,130:,:]
    u10=fv2.U10[hour,130:,:]
    psl=fv2.PSL[hour,130:,:]/100
    clat= fv_lin_centers_2.lat.values[np.where(np.logical_and(fv_lin_centers_2.label.values==label,fv_lin_centers_2.datetime.values==fv2.time[b].values))[0][0]]
    clon= fv_lin_centers_2.lon.values[np.where(np.logical_and(fv_lin_centers_2.label.values==label,fv_lin_centers_2.datetime.values==fv2.time[b].values))[0][0]]
    #print(clon,clat)
    if clat<70:
        continue
    files=[i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h2.*-{m}-{d}-00000.nc'.format(m=month,d=day))]
    files.sort()
    files1=[i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h1.*-{m}-{d}-00000.nc'.format(m=month,d=day))]
    files1.sort()
    temps=[]
    tpws=[]
    wind=[]
    psls=[]
    pvus=[]
    wind_850=[]
    for i,file in enumerate(files):
        fv2=xr.open_dataset(file)
        fv1=xr.open_dataset(files1[i])
        for b in [0,6,12,18]:
            if file=='/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h2.2013-01-01-00000.nc' and b==0:
                continue
            pressure=(fv1.hyam*fv1.P0+fv1.hybm*fv2.PS)/100
            ilev=fv1.ilev
            lat=fv1.lat.values[130:]
            lon=fv1.lon.values
            campressl=pressure.transpose("time", "lev", "lat", "lon") # hybrid coord press levels
            theta=weather_modules.temp_to_theta(fv1.T, campressl*100)
            theta_pres=hs2p(theta[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,ilev.data)
            u_pres=hs2p(fv1.U[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,ilev.data)
            v_pres=hs2p(fv1.V[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,ilev.data)
            t_pres=hs2p(fv1.T[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,ilev.data)
            q_pres=hs2p(fv1.Q[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,ilev.data)
            #z_pres=hs2p(fv3.Z3[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,ilev.data)
            press_levels3d = np.tile(ilev.data[:, np.newaxis, np.newaxis], (1, 62, 288))
            plev=np.array([300])
            #z_pres_300=hs2p(fv3.Z3[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,plev)/10
            epv=weather_modules.epv_sphere(theta_pres,press_levels3d*100,u_pres,v_pres,lat,lon)
            pvu_thet=weather_modules.interp2pv(epv, theta_pres, 2)
            plev=np.array([850])
            u_pres_850=hs2p(fv1.U[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,plev)
            v_pres_850=hs2p(fv1.V[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,plev)
            del(campressl)
            wind_850.append(np.sqrt(u_pres_850**2+v_pres_850**2))
            temps.append(fv2.T850[b,130:,:])
            tpws.append(fv2.TMQ[b,130:,:])
            wind.append(fv2.U10[b,130:,:])
            psls.append(fv2.PSL[b,130:,:])
            pvus.append(pvu_thet)
    temp_mean=np.mean(np.array(temps),axis=0)
    tpw_mean=np.mean(np.array(tpws),axis=0)
    u10_mean=np.mean(np.array(wind),axis=0)
    psl_mean=np.mean(np.array(psls),axis=0)
    u850_mean=np.mean(np.array(wind_850),axis=0)
    pvu_mean=np.mean(np.array(pvus),axis=0)
        
    psl_mean=np.mean((fv2.PSL)[:,130:,:],axis=0)/100
    slp_section=[]
    tpw_section=[]
    u10_section=[]
    psl_section=[]
    u850_section=[]
    pvu_section=[]
    center_section=[]
    non_center_section=[]
    for count, s in enumerate(slp):
        for c2, t in enumerate(s):
                dist = findDistance(clat,clon,t.lat,t.lon)
     #   print(dist)
                if dist <=1200:
                    slp_section.append(t-temp_mean[count,c2])
                    tpw_section.append(tpw[count,c2].values-tpw_mean[count,c2])
                    u10_section.append(u10[count,c2].values-u10_mean[count,c2])
                    psl_section.append(psl[count,c2].values-psl_mean[count,c2])
                    u850_section.append(u850[count,c2]-u850_mean[count,c2])
                    pvu_section.append(pvu_theta[count,c2]-pvu_mean[count,c2])
                    if dist <=500:
                       center_section.append(t)
                    else:
                        non_center_section.append(t)
    x = np.empty(len(slp_section))
    y = np.empty(len(slp_section))
    slp_v=np.empty(len(slp_section))
    ydiff = []
    xdiff = []
    for i in range(len(x)):
        slp_point = (slp_section[i].lat.values, slp_section[i].lon.values)
        center = (clat,clon)
        y_point = (slp_section[i].lat.values,clon )
        y[i] = distance.great_circle(center,y_point).km
        slp_v[i] = slp_section[i].values
        r = distance.great_circle(center,slp_point).km
        if y[i]>r:
            y[i]=r
        x[i] = np.sqrt(distance.great_circle(center,slp_point).km**2 - y[i]**2 )
        ydiff.append((slp_section[i].lat.values-clat) <= 0)
            #clon_180= ((clon + 180) % 360) - 180
            #slp_180= (( slp_section[i].lon.values+ 180) % 360) - 180
        xdiff.append((slp_section[i].lon.values-clon) <=0)
        if clon > 300 and slp_section[i].lon.values <60:
            xdiff[i]=False
        if clon < 60 and slp_section[i].lon.values >300:
            xdiff[i]=True
        if xdiff[i]:
            x[i] = -x[i]
        if ydiff[i]:
            y[i] = -y[i]
    start=(np.where(np.logical_and(fv_lin_centers_2.datetime.values==slp.time.values-timedelta(hours=1),fv_lin_centers_2.label.values==label))[0])
    end=(np.where(np.logical_and(fv_lin_centers_2.datetime.values==slp.time.values+timedelta(hours=6),fv_lin_centers_2.label.values==label))[0])
    if end: 
        #inds.append(a)
#print(inds)
        #zipped_lists = zip(xp, yp,slp_v)
        zipped_lists_n = zip(x, y,slp_v)
        zipped_lists_n1 = zip(x, y,np.array(tpw_section))
        zipped_lists_n2 = zip(x, y,np.array(u10_section))
        zipped_lists_n3 = zip(x, y,np.array(psl_section))
        zipped_lists_n4 = zip(x,y,np.array(u850_section))
        zipped_lists_n5 = zip(x,y,np.array(pvu_section))
            #sorted_pairs = sorted(zipped_lists)
        sorted_pairs_n = sorted(zipped_lists_n)
        sorted_pairs_n1 = sorted(zipped_lists_n1)
        sorted_pairs_n2 = sorted(zipped_lists_n2)
        sorted_pairs_n3 = sorted(zipped_lists_n3)
        sorted_pairs_n4 = sorted(zipped_lists_n4)
        sorted_pairs_n5 = sorted(zipped_lists_n5)

            #tuples = zip(*sorted_pairs)
            #list1, list2,list3 = [ list(tuple) for tuple in  tuples]
            
        tuples = zip(*sorted_pairs_n)
        tuples1 = zip(*sorted_pairs_n1)
        tuples2=zip(*sorted_pairs_n2)
        tuples3=zip(*sorted_pairs_n3)
        tuples4=zip(*sorted_pairs_n4)
        tuples5=zip(*sorted_pairs_n5)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        list4, list5,list6 = [ list(tuple) for tuple in  tuples1]
        list7, list8,list9 = [ list(tuple) for tuple in  tuples2]
        list10, list11,list12 = [ list(tuple) for tuple in  tuples3]
        list13, list14,list15 = [ list(tuple) for tuple in  tuples4]
        list16, list17,list18 = [ list(tuple) for tuple in  tuples5]
        ind_x=np.where(slp.lon==clon)[0][0]
        ind_y=np.where(slp.lat==clat)[0][0]
        if np.logical_and(np.logical_and(np.mean(np.array(center_section)) > np.mean(np.array(non_center_section)),slp[ind_y][ind_x]>0),np.mean(np.array(center_section))>0):
            warm_core.append([list1,list2,list3,list6,list9,list12,list15,list18])
        
        cold_core.append([list1,list2,list3,list6,list9,list12,list15,list18])
x = np.arange(-1200, 1201, 150)
y = np.arange(-1200, 1201, 150)
xv, yv = np.meshgrid(x, y)
list1=[row[0] for row in cold_core]
list2=[row[1] for row in cold_core]
list3=[row[2] for row in cold_core]
list4=[row[3] for row in cold_core]
list5=[row[4] for row in cold_core]
list6=[row[5] for row in cold_core]
list7=[row[6] for row in cold_core]
list8=[row[7] for row in cold_core]
print(len(list3))
list_1 = [item for row in list1 for item in row]
list_2 = [item for row in list2 for item in row]
list_3 = [item for row in list3 for item in row]
list_4 = [item for row in list4 for item in row]
list_5 = [item for row in list5 for item in row]
list_6 = [item for row in list6 for item in row]
list_7 = [item for row in list7 for item in row]
list_8 = [item for row in list8 for item in row]

        #means.append(np.nanmean(np.array(list_3)[ii[kk]]))
#print(len(list3))
#means=np.array(means).reshape(len(x)-1,-1)
#print(means)
#fig, ax = plt.subplots()
#cs=plt.pcolormesh(xv,yv,means)
#cs=plt.scatter(list_1,list_2,list_3)
#cbar = plt.colorbar(cs)
#plt.savefig('fv_lin_min_anom_150_all_test_north_2.png')
#znew0 = griddata((list_1,list_2),list_3, (xv, yv), method='linear',fill_value=np.nan)
#print(list1[0])
fig, ax = plt.subplots()
#plt.colorbar(cs)
#cs = plt.imshow(znew0, origin='lower',extent=[x[0], x[-1], y[0], y[-1]],cmap='bwr')
#cs=sns.kdeplot(np.array(list_3),x=x,y=y,weights=np.array(list_3),cmap='bwr')
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_3))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='bwr',vmin=-10,vmax=10,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='[K]',extend='both')
#plt.savefig('fv_25_max3672_150_cold_test_north_temp_70.png')

fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_4))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='BrBG',vmin=-5,vmax=5,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label=r'[kg/m$^2$]',extend='both')
#plt.savefig('fv_25_max3672_150_cold_test_north_tpw_70.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_5))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-5,vmax=5,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='m/s',extend='both')
#plt.savefig('fv_25_max3672_150_cold_test_north_u10_70.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_6))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-40,vmax=40,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='hPa',extend='both')
#plt.savefig('fv_25_max3672_150_cold_test_north_psl_70.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_7))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-5,vmax=5,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='m/s',extend='both')
#plt.savefig('fv_25_min_150_cold_test_north_u850.png')
fig, ax = plt.subplots()
nan_mask=np.isnan(np.array(list_8))
H, xedges, yedges = np.histogram2d(np.array(list_1)[~nan_mask],np.array(list_2)[~nan_mask], bins = [x, y], weights = np.array(list_8)[~nan_mask])
H_counts, xedges, yedges = np.histogram2d(np.array(list_1)[~nan_mask],np.array(list_2)[~nan_mask], bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-20,vmax=20,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='K',extend='both')
plt.savefig('fv_25_max3562_150_cold_test_north_pvu_70.png')
'''fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_6))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-20,vmax=20,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='m/s',extend='both')
plt.savefig('fv_25_min_anom_150_all_test_north_12.png')'''

#H_counts, xedges, yedges = np.histogram2d(x=np.array(list_1), y=np.array(list_2), bins = [16, 16]) 
#cs=plt.hist2d(x=np.array(list_1), y=np.array(list_2),bins=[16, 16],range=[[-1200, 1200], [-1200, 1200]],weights=np.array(list_3),cmap='bwr')
#cs=plt.tripcolor(np.array(list_1),np.array(list_2),np.array(list_3),cmap='jet')
#cbar = plt.colorbar(cs[3])
#plt.savefig('fv_lin_min_anom_150_all_test_north_2.png')
list1=[row[0] for row in warm_core]
list2=[row[1] for row in warm_core]
list3=[row[2] for row in warm_core]
list4=[row[3] for row in warm_core]
list5=[row[4] for row in warm_core]
list6=[row[5] for row in warm_core]
list7=[row[6] for row in warm_core]
list8=[row[7] for row in warm_core]
print(len(list3))
'''list_1 = [item for row in list1 for item in row]
list_2 = [item for row in list2 for item in row]
list_3 = [item for row in list3 for item in row]
list_4 = [item for row in list4 for item in row]
list_5 = [item for row in list5 for item in row]
list_6 = [item for row in list6 for item in row]
list_7 = [item for row in list7 for item in row]
list_8 = [item for row in list8 for item in row]
fig, ax = plt.subplots()
#plt.colorbar(cs)
#cs = plt.imshow(znew0, origin='lower',extent=[x[0], x[-1], y[0], y[-1]],cmap='bwr')
#cs=sns.kdeplot(np.array(list_3),x=x,y=y,weights=np.array(list_3),cmap='bwr')

H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_3))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='bwr',vmin=-5,vmax=5,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='[K]',extend='both')
plt.savefig('fv_25_max_150_warm_test_north_temp_70.png')

fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_4))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='BrBG',vmin=-3,vmax=3,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label=r'[kg/m$^2$]',extend='both')
plt.savefig('fv_25_max_150_warm_test_north_tpw_70.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_5))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-5,vmax=5,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='m/s',extend='both')
plt.savefig('fv_25_max_150_warm_test_north_u10_70.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_6))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-20,vmax=20,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='hPa',extend='both')
plt.savefig('fv_25_max_150_warm_test_north_psl_70.png')
fig, ax = plt.subplots()
nan_mask=np.isnan(np.array(list_8))
H, xedges, yedges = np.histogram2d(np.array(list_1)[~nan_mask],np.array(list_2)[~nan_mask], bins = [x, y], weights = np.array(list_8)[~nan_mask])
H_counts, xedges, yedges = np.histogram2d(np.array(list_1)[~nan_mask],np.array(list_2)[~nan_mask], bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-20,vmax=20,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='K',extend='both')
plt.savefig('fv_25_max_150_warm_test_north_pvu_70.png')'''