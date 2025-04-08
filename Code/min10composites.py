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
from scipy.interpolate import griddata
import nc_time_axis
import os
os.environ['ESMFMKFILE']='/glade/work/ajacobs/conda-envs/asp2023/lib/esmf.mk'
import xesmf as xe
import pandas as pd
import seaborn as sns
import warnings
from RegriddingFunctions import hs2p
import weather_modules
warnings.filterwarnings("ignore")
cyc=xr.open_dataset('all_nudge_abl_cyclone_centers.nc')
end=int(cyc.label.values[-1][-4:])
cyc=cyc.sortby(['label','slp'])
warm_core=[]
cold_core=[]
spring=[]
summer=[]
winter=[]
fall=[]
def findDistance(lat1,lon1,lat2,lon2):
        center = (lat1,lon1)
        point = (lat2,lon2)
        dist = distance.great_circle(center,point).km
        return dist

for a in [3562]: #3566,3562,288,642,1666,3672#3566,361,4220,4709,4716,4724,843,845
    '''2013-07-24 06:00:00 cyclone_288
2014-05-16 00:00:00 cyclone_642
2014-09-08 12:00:00 cyclone_806
2015-09-30 12:00:00 cyclone_1345
2016-06-08 12:00:00 cyclone_1666
2016-08-15 06:00:00 cyclone_1747
2018-06-07 12:00:00 cyclone_2655
2018-07-07 12:00:00 cyclone_2697
2018-07-19 00:00:00 cyclone_2708
2020-05-09 06:00:00 cyclone_3562
2020-05-12 18:00:00 cyclone_3566
2020-07-28 06:00:00 cyclone_3672
2020-08-14 12:00:00 cyclone_3699'''
    label = 'cyclone_'+str(a).zfill(3)
    time= cyc.datetime.values[np.where(cyc.label.values==label)[0][0]]
    en=(np.where(np.logical_and(cyc.datetime.values==time+timedelta(hours=6),cyc.label.values==label))[0])
    #time=time+timedelta(hours=6)
    print(en)
    time+=timedelta(hours=6)
    #clat= cyc.lat.values[np.where(cyc.label.values==label)[0][0]]
    #if not en:
    #    continue
    clat=cyc.lat.values[en]
    clon=cyc.lon.values[en]
    cslp=cyc.slp.values[en]
    #print(clon)
    #clon= cyc.lon.values[np.where(cyc.label.values==label)[0][0]]
    #cslp=cyc.slp.values[np.where(cyc.label.values==label)[0][0]]
    year=str(time.year)
    month=str(time.month).zfill(2)
    day=str(time.day).zfill(2)
    files=[i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h2.*-{m}-{d}-00000.nc'.format(m=month,d=day))]
    files.sort()
    files1=[i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h1.*-{m}-{d}-00000.nc'.format(m=month,d=day))]
    files1.sort()
    temps=[]
    tpws=[]
    wind=[]
    psls=[]
    wind_850=[]
    pvus=[]
    pvqs=[]
    pvwinds=[]
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
            plev=np.array([850])
            u_pres_850=hs2p(fv1.U[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,plev)
            v_pres_850=hs2p(fv1.V[b:b+1,:,130:,:].data,campressl[b:b+1,:,130:,:].data,plev)
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
            del(campressl)
            epv=weather_modules.epv_sphere(theta_pres,press_levels3d*100,u_pres,v_pres,lat,lon)
            pvu_theta=weather_modules.interp2pv(epv, theta_pres, 2)
            pvu_q=weather_modules.interp2pv(epv, q_pres, 2)
            pvu_u=weather_modules.interp2pv(epv, u_pres, 2)
            pvu_v=weather_modules.interp2pv(epv, v_pres, 2)
            wind_850.append(np.sqrt(u_pres_850**2+v_pres_850**2))
            temps.append(fv2.T850[b,130:,:])
            tpws.append(fv2.TMQ[b,130:,:])
            wind.append(fv2.U10[b,130:,:])
            psls.append(fv2.PSL[b,130:,:])
            pvus.append(pvu_theta)
            pvqs.append(pvu_q*1000)
            pvwinds.append(np.sqrt(pvu_u**2+pvu_v**2))
    temp_mean=np.mean(np.array(temps),axis=0)
    tpw_mean=np.mean(np.array(tpws),axis=0)
    u10_mean=np.mean(np.array(wind),axis=0)
    psl_mean=np.mean(np.array(psls),axis=0)
    u850_mean=np.mean(np.array(wind_850),axis=0)
    pvu_mean=np.mean(np.array(pvus),axis=0)
    pvq_mean=np.mean(np.array(pvqs),axis=0)
    pvwinds_mean=np.mean(np.array(pvwinds),axis=0)
    file='/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h2.{y}-{m}-{d}-00000.nc'.format(y=year,m=month,d=day)
    fv2=xr.open_dataset(file)
    file1='/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h1.{y}-{m}-{d}-00000.nc'.format(y=year,m=month,d=day)
    fv1=xr.open_dataset(file1)
    hour=int(time.hour)
    if file=='/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h2.2013-01-01-00000.nc' and hour==0:
        continue
    pressure=(fv1.hyam*fv1.P0+fv1.hybm*fv2.PS)/100
    ilev=fv1.ilev
    lat=fv1.lat.values[130:]
    lon=fv1.lon.values
    campressl=pressure.transpose("time", "lev", "lat", "lon") # hybrid coord press levels
    plev=np.array([850])
    u_pres_850=hs2p(fv1.U[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,plev)
    v_pres_850=hs2p(fv1.V[hour:hour+1,:,130:,:].data,campressl[hour:hour+1,:,130:,:].data,plev)
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
    pvu_q=weather_modules.interp2pv(epv, q_pres, 2)*1000
    pvu_u=weather_modules.interp2pv(epv, u_pres, 2)
    pvu_v=weather_modules.interp2pv(epv, v_pres, 2)
    del(campressl)
    u850=np.sqrt(u_pres_850**2+v_pres_850**2)
    pvwind=np.sqrt(pvu_u**2+pvu_v**2)
    t=fv2.T850[hour,130:,:]
    slp=(fv2.T850)[hour,130:,:]
    tpw=fv2.TMQ[hour,130:,:]
    u10=fv2.U10[hour,130:,:]
    psl=fv2.PSL[hour,130:,:]/100
    slp_section=[]
    tpw_section=[]
    u10_section=[]
    psl_section=[]
    u850_section=[]
    pvu_section=[]
    pvq_section=[]
    pvwind_section=[]
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
                pvq_section.append(pvu_q[count,c2]-pvq_mean[count,c2])
                pvwind_section.append(pvwind[count,c2]-pvwinds_mean[count,c2])
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
    start=(np.where(np.logical_and(cyc.datetime.values==slp.time.values-timedelta(hours=6),cyc.label.values==label))[0])
    end=(np.where(np.logical_and(cyc.datetime.values==slp.time.values+timedelta(hours=6),cyc.label.values==label))[0])
    if True:
        if not start:
            start_lat = clat
            start_lon = clon
        start_lat = cyc.lat.values[start]
        start_lon = cyc.lon.values[start]
        end_lat = cyc.lat.values[end]
        end_lon = cyc.lon.values[end]
        
        start = (start_lat,start_lon)
        end = (end_lat,end_lon)
        y_point = (end_lat,start_lon)
        y_center = distance.great_circle(start,y_point).km 
        x_center = np.sqrt(distance.great_circle(start,end).km**2 - y_center**2)
# NW quadrant
        if (end_lon - start_lon < 0) and (end_lat - start_lat >= 0):
            heading = np.pi - np.arctan( y_center / x_center )
    # NE quadrant
        elif (end_lon - start_lon >= 0) and (end_lat - start_lat >= 0):
            heading = np.arctan2( y_center, x_center )
    # SE quadrant
        elif (end_lon - start_lon >= 0) and (end_lat - start_lat < 0):
            heading = -np.arctan( y_center / x_center )
    # SW quadrant
        else:
            heading = -np.pi + np.arctan(y_center / x_center)
        xp = x * np.cos(heading-np.pi/2) + y * np.sin(heading-np.pi/2)
        yp = -x * np.sin(heading-np.pi/2) + y * np.cos(heading-np.pi/2)
        zipped_lists = zip(xp, yp,slp_v)
        zipped_lists_n = zip(x, y,slp_v)
        zipped_lists_n1 = zip(x, y,np.array(tpw_section))
        zipped_lists_n2 = zip(x, y,np.array(u10_section))
        zipped_lists_n3 = zip(x, y,np.array(psl_section))
        zipped_lists_n4 = zip(x,y,np.array(u850_section))
        zipped_lists_n5 = zip(x,y,np.array(pvu_section))
        zipped_lists_n6 = zip(x,y,np.array(pvq_section))
        zipped_lists_n7 = zip(x,y,np.array(pvwind_section))
            #sorted_pairs = sorted(zipped_lists)
        sorted_pairs_n = sorted(zipped_lists_n)
        sorted_pairs_n1 = sorted(zipped_lists_n1)
        sorted_pairs_n2 = sorted(zipped_lists_n2)
        sorted_pairs_n3 = sorted(zipped_lists_n3)
        sorted_pairs_n4 = sorted(zipped_lists_n4)
        sorted_pairs_n5 = sorted(zipped_lists_n5)
        sorted_pairs_n6 = sorted(zipped_lists_n6)
        sorted_pairs_n7 = sorted(zipped_lists_n7)

            #tuples = zip(*sorted_pairs)
            #list1, list2,list3 = [ list(tuple) for tuple in  tuples]
            
        tuples = zip(*sorted_pairs_n)
        tuples1 = zip(*sorted_pairs_n1)
        tuples2=zip(*sorted_pairs_n2)
        tuples3=zip(*sorted_pairs_n3)
        tuples4=zip(*sorted_pairs_n4)
        tuples5=zip(*sorted_pairs_n5)
        tuples6=zip(*sorted_pairs_n6)
        tuples7=zip(*sorted_pairs_n7)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        list4, list5,list6 = [ list(tuple) for tuple in  tuples1]
        list7, list8,list9 = [ list(tuple) for tuple in  tuples2]
        list10, list11,list12 = [ list(tuple) for tuple in  tuples3]
        list13, list14,list15 = [ list(tuple) for tuple in  tuples4]
        list16, list17,list18 = [ list(tuple) for tuple in  tuples5]
        list19, list20,list21 = [ list(tuple) for tuple in  tuples6]
        list22, list23,list24 = [ list(tuple) for tuple in  tuples7]
        print(slp.lon)
        ind_x=np.where(slp.lon==clon)[0][0]
        ind_y=np.where(slp.lat==clat)[0][0]
        if np.logical_and(np.mean(np.array(center_section)) > np.mean(np.array(non_center_section)),(slp[ind_y][ind_x]-temp_mean[ind_y][ind_x])>0):
            warm_core.append([list1,list2,list3,list6,list9,list12,list15,list18])
            #else:
        cold_core.append([list1,list2,list3,list6,list9,list12,list15,list18,list21,list24])
        if int(month) in [5,6,7,8,9]:
            summer.append([list1,list2,list3,list6,list9,list12,list15,list18,list21,list24])
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
list9=[row[8] for row in cold_core]
list10=[row[9] for row in cold_core]
print(len(list3))
list_1 = [item for row in list1 for item in row]
list_2 = [item for row in list2 for item in row]
list_3 = [item for row in list3 for item in row]
list_4 = [item for row in list4 for item in row]
list_5 = [item for row in list5 for item in row]
list_6 = [item for row in list6 for item in row]
list_7 = [item for row in list7 for item in row]
list_8 = [item for row in list8 for item in row]
list_9 = [item for row in list9 for item in row]
list_10 = [item for row in list10 for item in row]

        #means.append(np.nanmean(np.array(list_3)[ii[kk]]))
print(len(list3))
#means=np.array(means).reshape(len(x)-1,-1)
#print(means)
#fig, ax = plt.subplots()
#cs=plt.pcolormesh(xv,yv,means)
#cs=plt.scatter(list_1,list_2,list_3)
#cbar = plt.colorbar(cs)
#plt.savefig('fv_lin_min_anom_150_all_test_north_2.png')
#znew0 = griddata((list_1,list_2),list_3, (xv, yv), method='linear',fill_value=np.nan)
#print(list1[0])
'''fig, ax = plt.subplots()
#plt.colorbar(cs)
#cs = plt.imshow(znew0, origin='lower',extent=[x[0], x[-1], y[0], y[-1]],cmap='bwr')
#cs=sns.kdeplot(np.array(list_3),x=x,y=y,weights=np.array(list_3),cmap='bwr')
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_3))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='bwr',vmin=-10,vmax=10,#vmin=250,vmax=290,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='[K]',extend='both')
plt.savefig('fv_25_minplus6h_150_min3562_test_north_temp.png')

fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_4))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='BrBG',vmin=-5,vmax=5,#cmap='viridis',vmin=0,vmax=50,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label=r'[kg/m$^2$]',extend='both')
plt.savefig('fv_25_minplus6h_150_min3562_test_north_tpw.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_5))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower', cmap='PRGn',vmin=-10,vmax=10, #cmap='YlOrRd',vmin=0,vmax=25,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='m/s',extend='both')
#plt.savefig('fv_25_minplus1_150_min3562_test_north_u10.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_6))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower', cmap='PRGn',vmin=-40,vmax=40, #cmap='plasma',vmin=960,vmax=1010,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='hPa',extend='both')
plt.savefig('fv_25_minplus6h_150_min3562_test_north_psl.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_7))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-5,vmax=5,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='m/s',extend='both')
#plt.savefig('fv_25_min_anom_150_all_test_north_u850.png')
fig, ax = plt.subplots()
nan_mask=np.isnan(np.array(list_8))
H, xedges, yedges = np.histogram2d(np.array(list_1)[~nan_mask],np.array(list_2)[~nan_mask], bins = [x, y], weights = np.array(list_8)[~nan_mask])
H_counts, xedges, yedges = np.histogram2d(np.array(list_1)[~nan_mask],np.array(list_2)[~nan_mask], bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-20,vmax=20,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='K',extend='both')
plt.savefig('fv_25_minplus6h_150_min5_test_north_pvu.png')
fig, ax = plt.subplots()
nan_mask=np.isnan(np.array(list_9))
H, xedges, yedges = np.histogram2d(np.array(list_1)[~nan_mask],np.array(list_2)[~nan_mask], bins = [x, y], weights = np.array(list_9)[~nan_mask])
H_counts, xedges, yedges = np.histogram2d(np.array(list_1)[~nan_mask],np.array(list_2)[~nan_mask], bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='BrBG',vmin=-3,vmax=3,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='m/s',extend='both')
#plt.savefig('fv_25_min_150_min5_test_north_pvq.png')
fig, ax = plt.subplots()
nan_mask=np.isnan(np.array(list_10))
H, xedges, yedges = np.histogram2d(np.array(list_1)[~nan_mask],np.array(list_2)[~nan_mask], bins = [x, y], weights = np.array(list_10)[~nan_mask])
H_counts, xedges, yedges = np.histogram2d(np.array(list_1)[~nan_mask],np.array(list_2)[~nan_mask], bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-5,vmax=5,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='m/s',extend='both')
#plt.savefig('fv_25_min_150_min5_test_north_upv.png')
fig, ax = plt.subplots()
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
print(len(list3))
'''list_1 = [item for row in list1 for item in row]
list_2 = [item for row in list2 for item in row]
list_3 = [item for row in list3 for item in row]
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
#plt.savefig('fv_25_min_anom_150_warm_test_north_6.png')
#print(list3[0])
#znew0 = griddata((list_1,list_2),list_3, (xv, yv), method='linear')
#fig, ax = plt.subplots()

#cs = plt.imshow(znew0, origin='lower',extent=[x[0], x[-1], y[0], y[-1]],cmap='bwr')

#cbar = plt.colorbar(cs)
#plt.savefig('fv_lin_min_anom_200_all_warm_north_2.png')

x = np.arange(-1200, 1201, 150)
y = np.arange(-1200, 1201, 150)
xv, yv = np.meshgrid(x, y)
list1=[row[0] for row in summer]
list2=[row[1] for row in summer]
list3=[row[2] for row in summer]
list4=[row[3] for row in summer]
list5=[row[4] for row in summer]
list6=[row[5] for row in summer]
list7=[row[6] for row in summer]
print(len(list3))
list_1 = [item for row in list1 for item in row]
list_2 = [item for row in list2 for item in row]
list_3 = [item for row in list3 for item in row]
list_4 = [item for row in list4 for item in row]
list_5 = [item for row in list5 for item in row]
list_6 = [item for row in list6 for item in row]
list_7 = [item for row in list7 for item in row]

        #means.append(np.nanmean(np.array(list_3)[ii[kk]]))
print(len(list3))
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

plt.imshow(H.T, origin='lower',  cmap='bwr',vmin=-5,vmax=5,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='[K]',extend='both')
plt.savefig('fv_25_min_anom_150_summer_test_north_temp.png')

fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_4))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='BrBG',vmin=-3,vmax=3,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label=r'[kg/m$^2$]',extend='both')
plt.savefig('fv_25_min_anom_150_summer_test_north_tpw.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_5))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-5,vmax=5,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='m/s',extend='both')
plt.savefig('fv_25_min_anom_150_summer_test_north_u10.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_6))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-20,vmax=20,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='hPa',extend='both')
plt.savefig('fv_25_min_anom_150_summer_test_north_psl.png')
fig, ax = plt.subplots()
H, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y], weights = np.array(list_7))
H_counts, xedges, yedges = np.histogram2d(np.array(list_1),np.array(list_2), bins = [x, y])
H = H/H_counts

plt.imshow(H.T, origin='lower',  cmap='PRGn',vmin=-5,vmax=5,
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
plt.colorbar(label='m/s',extend='both')
plt.savefig('fv_25_min_anom_150_summer_test_north_u850.png')'''
