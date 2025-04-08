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
from metpy.interpolate import cross_section
from RegriddingFunctions import hs2p


files = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h1.2022-07-0*.nc')]
files.sort()
files=files[6:7]#6:7
files2 = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h2.2022-07-0*.nc')]
files2.sort()
files2=files2[6:7]
files3 = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h3.2022-07-0*.nc')]
files3.sort()
files3=files3[6:7]
fv_25_centers=xr.open_dataset('2022_nudge_abl_cyclone_centers.nc')
files_a500 = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022.cam.h1.2022-07-0*.nc')]
files_a500.sort()
files_a500=files_a500[6:7]
files2_a500 = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022.cam.h2.2022-07-0*.nc')]
files2_a500.sort()
files2_a500=files2_a500[6:7]
files3_a500 = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022.cam.h3.2022-07-0*.nc')]
files3_a500.sort()
files3_a500=files3_a500[6:7]
files_at = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022.cam.h1.2022-07-0*.nc')]
files_at.sort()
files_at=files_at[6:7]
files2_at = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022.cam.h2.2022-07-0*.nc')]
files2_at.sort()
files2_at=files2_at[6:7]
files3_at = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022.cam.h3.2022-07-0*.nc')]
files3_at.sort()
files3_at=files3_at[6:7]
fv_17_centers=xr.open_dataset('2022_nudge_at_cyclone_centers.nc')
fv_20_centers=xr.open_dataset('2022_nudge_a500_cyclone_centers.nc')
era5_fv=xr.open_dataset('2022_era5_cyclone_centers.nc')
era5pres=xr.open_dataset('/glade/derecho/scratch/ajacobs/cyclone_nc_files/era5_fv_cyc2_pres.nc')
era5sfc=xr.open_dataset('/glade/derecho/scratch/ajacobs/cyclone_nc_files/era5_fv_cyc2_sfc.nc')

def findDistance(lat1,lon1,lat2,lon2):
        center = (lat1,lon1)
        point = (lat2,lon2)
        dist = distance.great_circle(center,point).km
        return dist
for k,file in enumerate(files):
    fv3=xr.open_dataset(files3[k])
    fv1=xr.open_dataset(file)
    fv2=xr.open_dataset(files2[k])
    fv3_at=xr.open_dataset(files3_at[k])
    fv1_at=xr.open_dataset(files_at[k])
    fv2_at=xr.open_dataset(files2_at[k])
    fv3_a500=xr.open_dataset(files3_a500[k])
    fv1_a500=xr.open_dataset(files_a500[k])
    fv2_a500=xr.open_dataset(files2_a500[k])
    for a in [6]: #18
        snow_flux =  (fv2.PRECSC+fv2.PRECSL*997*3.337e5)
    #(fv2.SHFLX+fv2.LHFLX+fv2.FLNS-fv2.FSNS+snow_flux)
    #slp=(fv2.U10)[a,130:,:]
        slp=fv2.PSL[a,130:,]/100
        tpw=fv2.TMQ[a,130:,]
        t850=fv2.T850[a,130:,]
        u10=fv2.U10[a,130:,:]
        pressure=(fv1.hyam*fv1.P0+fv1.hybm*fv2.PS[a:a+1])/100
        plev=np.array([500])
        campressl=pressure.transpose("time", "lev", "lat", "lon")
        omega=fv3.OMEGA[a:a+1,:,130:,:]*36
        o_pres=hs2p(omega[:,:,:,:].data,campressl[:,:,130:,:].data,plev)
        z_pres=hs2p(fv3.Z3[a:a+1,:,130:,:].data,campressl[:,:,130:,:].data,plev)/10
        del(campressl)
    
    
    
        clon=fv_25_centers.lon.values[np.where(np.logical_and(fv_25_centers.datetime.values==slp.time.values,fv_25_centers.label.values=='cyclone_233'))[0]]#cyclone_020
        clat=fv_25_centers.lat.values[np.where(np.logical_and(fv_25_centers.datetime.values==slp.time.values,fv_25_centers.label.values=='cyclone_233'))[0]]
        slp_section=[]
        tpw_section=[]
        t850_section=[]
        u10_section=[]
        o_section=[]
        z_section=[]
        indices=[]
        for count, s in enumerate(slp):
            for c2, t in enumerate(s):
                dist = findDistance(clat,clon,t.lat,t.lon)
     #   print(dist)
                if dist <=1200:
                    slp_section.append(t)
                    indices.append([count+130,c2])
                    tpw_section.append(tpw[count,c2].values)
                    t850_section.append(t850.values[count,c2])
                    u10_section.append(u10[count,c2].values)
                    o_section.append(o_pres[count,c2])
                    z_section.append(z_pres[count,c2])
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
            xdiff.append((slp_section[i].lon.values-clon) <=0)
            if clon > 300 and slp_section[i].lon.values <60:
                xdiff[i]=False
            if clon < 60 and slp_section[i].lon.values >300:
                xdiff[i]=True
            if xdiff[i]:
                x[i] = -x[i]
            if ydiff[i]:
                y[i] = -y[i]
        start_lat = fv_25_centers.lat.values[(np.where(np.logical_and(fv_25_centers.datetime.values==slp.time.values-timedelta(hours=1),fv_25_centers.label.values=='cyclone_233'))[0])]
        start_lon = fv_25_centers.lon.values[(np.where(np.logical_and(fv_25_centers.datetime.values==slp.time.values-timedelta(hours=1),fv_25_centers.label.values=='cyclone_233'))[0])]
        end_lat = fv_25_centers.lat.values[(np.where(np.logical_and(fv_25_centers.datetime.values==slp.time.values+timedelta(hours=1),fv_25_centers.label.values=='cyclone_233'))[0])]
        end_lon = fv_25_centers.lon.values[(np.where(np.logical_and(fv_25_centers.datetime.values==slp.time.values+timedelta(hours=1),fv_25_centers.label.values=='cyclone_233'))[0])]
        slp_at=fv2_at.PSL[a-5,130:,]/100
        tpw_at=fv2_at.TMQ[a-5,130:,]
        t850_at=fv2_at.T850[a-5,130:,]
        u10_at=fv2_at.U10[a-5,130:,:]
        pressure=(fv1_at.hyam*fv1_at.P0+fv1_at.hybm*fv2_at.PS[a-5:a-4])/100
        plev=np.array([500])
        campressl=pressure.transpose("time", "lev", "lat", "lon")
        omega_at=fv3_at.OMEGA[a-5:a-4,:,130:,:]*36
        o_pres_at=hs2p(omega_at[:,:,:,:].data,campressl[:,:,130:,:].data,plev)
        z_pres_at=hs2p(fv3_at.Z3[a-5:a-4,:,130:,:].data,campressl[:,:,130:,:].data,plev)/10
        del(campressl)
    
    
    
        clon_at=fv_17_centers.lon.values[np.where(np.logical_and(fv_17_centers.datetime.values==slp_at.time.values,fv_17_centers.label.values=='cyclone_247'))[0]]
        clat_at=fv_17_centers.lat.values[np.where(np.logical_and(fv_17_centers.datetime.values==slp_at.time.values,fv_17_centers.label.values=='cyclone_247'))[0]]#247
        slp_section_at=[]
        tpw_section_at=[]
        t850_section_at=[]
        u10_section_at=[]
        o_section_at=[]
        z_section_at=[]
        indices=[]
        for count, s in enumerate(slp_at):
            for c2, t in enumerate(s):
                dist = findDistance(clat_at,clon_at,t.lat,t.lon)
     #   print(dist)
                if dist <=1200:
                    slp_section_at.append(t)
                    indices.append([count+130,c2])
                    tpw_section_at.append(tpw_at[count,c2].values)
                    t850_section_at.append(t850_at.values[count,c2])
                    u10_section_at.append(u10_at[count,c2].values)
                    o_section_at.append(o_pres_at[count,c2])
                    z_section_at.append(z_pres_at[count,c2])
        x_at = np.empty(len(slp_section_at))
        y_at = np.empty(len(slp_section_at))
        slp_v_at=np.empty(len(slp_section_at))
        ydiff = []
        xdiff = []
        for i in range(len(x_at)):
            slp_point = (slp_section_at[i].lat.values, slp_section_at[i].lon.values)
            center = (clat_at,clon_at)
            y_point = (slp_section_at[i].lat.values,clon_at )
            y_at[i] = distance.great_circle(center,y_point).km
            slp_v_at[i] = slp_section_at[i].values
            r = distance.great_circle(center,slp_point).km
            x_at[i] = np.sqrt(distance.great_circle(center,slp_point).km**2 - y_at[i]**2 )
            ydiff.append((slp_section_at[i].lat.values-clat_at) <= 0)
            xdiff.append((slp_section_at[i].lon.values-clon_at) <=0)
            if clon_at > 300 and slp_section_at[i].lon.values <60:
                xdiff[i]=False
            if clon_at < 60 and slp_section_at[i].lon.values >300:
                xdiff[i]=True
            if xdiff[i]:
                x_at[i] = -x_at[i]
            if ydiff[i]:
                y_at[i] = -y_at[i]
        slp_a500=fv2_a500.PSL[a-6,130:,]/100
        tpw_a500=fv2_a500.TMQ[a-6,130:,]
        t850_a500=fv2_a500.T850[a-6,130:,]
        u10_a500=fv2_a500.U10[a-6,130:,:]
        pressure=(fv1_a500.hyam*fv1_a500.P0+fv1_a500.hybm*fv2_a500.PS[a-6:a-5])/100
        plev=np.array([500])
        campressl=pressure.transpose("time", "lev", "lat", "lon")
        omega_a500=fv3_a500.OMEGA[a-6:a-5,:,130:,:]*36
        o_pres_a500=hs2p(omega_a500[:,:,:,:].data,campressl[:,:,130:,:].data,plev)
        z_pres_a500=hs2p(fv3_a500.Z3[a-6:a-5,:,130:,:].data,campressl[:,:,130:,:].data,plev)/10
        del(campressl)
    
    
    
        clon_a500=fv_20_centers.lon.values[np.where(np.logical_and(fv_20_centers.datetime.values==slp_a500.time.values,fv_20_centers.label.values=='cyclone_226'))[0]]
        clat_a500=fv_20_centers.lat.values[np.where(np.logical_and(fv_20_centers.datetime.values==slp_a500.time.values,fv_20_centers.label.values=='cyclone_226'))[0]]
        slp_section_a500=[]
        tpw_section_a500=[]
        t850_section_a500=[]
        u10_section_a500=[]
        o_section_a500=[]
        z_section_a500=[]
        indices=[]
        for count, s in enumerate(slp_a500):
            for c2, t in enumerate(s):
                dist = findDistance(clat_a500,clon_a500,t.lat,t.lon)
     #   print(dist)
                if dist <=1200:
                    slp_section_a500.append(t)
                    indices.append([count+130,c2])
                    tpw_section_a500.append(tpw_a500[count,c2].values)
                    t850_section_a500.append(t850_a500.values[count,c2])
                    u10_section_a500.append(u10_a500[count,c2].values)
                    o_section_a500.append(o_pres_a500[count,c2])
                    z_section_a500.append(z_pres_a500[count,c2])
        x_a500 = np.empty(len(slp_section_a500))
        y_a500 = np.empty(len(slp_section_a500))
        slp_v_a500=np.empty(len(slp_section_a500))
        ydiff = []
        xdiff = []
        for i in range(len(x_a500)):
            slp_point = (slp_section_a500[i].lat.values, slp_section_a500[i].lon.values)
            center = (clat_a500,clon_a500)
            y_point = (slp_section_a500[i].lat.values,clon_a500 )
            y_a500[i] = distance.great_circle(center,y_point).km
            slp_v_a500[i] = slp_section_a500[i].values
            r = distance.great_circle(center,slp_point).km
            x_a500[i] = np.sqrt(distance.great_circle(center,slp_point).km**2 - y_a500[i]**2 )
            ydiff.append((slp_section_a500[i].lat.values-clat_a500) <= 0)
            xdiff.append((slp_section_a500[i].lon.values-clon_a500) <=0)
            if clon_a500 > 300 and slp_section_a500[i].lon.values <60:
                xdiff[i]=False
            if clon_a500 < 60 and slp_section_a500[i].lon.values >300:
                xdiff[i]=True
            if xdiff[i]:
                x_a500[i] = -x_a500[i]
            if ydiff[i]:
                y_a500[i] = -y_a500[i]
        slp_era5=era5sfc.msl[34,130:,:]/100 #34 #14
        tpw_era5=era5sfc.tcwv[34,130:,:]
        t850_era5=era5pres.t[26,-7,130:,:] #26 #6
        o_era5=era5pres.w[26,-16,130:,:]*36
        u10_era5=np.sqrt(era5sfc.u10[34,130:,:]**2+era5sfc.v10[34,130:,:]**2)
        lat=era5pres.lat.values[130:]
        lon=era5pres.lon.values
        clon_era5=era5_fv.lon.values[np.where(np.logical_and(era5_fv.datetime.values==slp_era5.time.values,era5_fv.label.values=='cyclone_241'))[0]] #cyclone_159
        clat_era5=era5_fv.lat.values[np.where(np.logical_and(era5_fv.datetime.values==slp_era5.time.values,era5_fv.label.values=='cyclone_241'))[0]]
        slp_section_era5=[]
        tpw_section_era5=[]
        t850_section_era5=[]
        o_section_era5=[]
        u10_section_era5=[]
        indices=[]
        for count, s in enumerate(slp_era5):
            for c2, t in enumerate(s):
                dist = findDistance(clat_era5,clon_era5,lat[count],lon[c2])
     #   print(dist)
                if dist <=1200:
                    slp_section_era5.append(t)
                    indices.append([count,c2])
                    tpw_section_era5.append(tpw_era5[count,c2].values)
                    t850_section_era5.append(t850_era5.values[count,c2])
                    o_section_era5.append(o_era5.values[count,c2])
                    u10_section_era5.append(u10_era5[count,c2])
        x_era5 = np.empty(len(slp_section_era5))
        y_era5 = np.empty(len(slp_section_era5))
        slp_v_era5=np.empty(len(slp_section_era5))
        ydiff = []
        xdiff = []
        for i in range(len(x_era5)):
            slp_point = (lat[indices[i][0]], lon[indices[i][1]])
            center = (clat_era5,clon_era5)
            y_point = (lat[indices[i][0]],clon_era5 )
            y_era5[i] = distance.great_circle(center,y_point).km
            slp_v_era5[i] = slp_section_era5[i].values
            r = distance.great_circle(center,slp_point).km
            x_era5[i] = np.sqrt(distance.great_circle(center,slp_point).km**2 - y_era5[i]**2 )
            ydiff.append((lat[indices[i][0]]-clat_era5) <= 0)
            xdiff.append((lon[indices[i][1]]-clon_era5) <=0)
            if clon_era5 > 300 and lon[indices[i][1]] <60:
                xdiff[i]=False
            if clon_era5 < 60 and lon[indices[i][1]] >300:
                xdiff[i]=True
            if xdiff[i]:
                x_era5[i] = -x_era5[i]
            if ydiff[i]:
                y_era5[i] = -y_era5[i]
        
        zipped_lists = zip(x_era5, y_era5,np.array(t850_section_era5))
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear')
        fig = plt.figure(figsize=[21,14])
        gs = gridspec.GridSpec(2,5, wspace= 0.5,hspace=None,width_ratios=[1,1,1,1,0.05])
        title='850 hPa Temperature and Total Precipatble Water'+'\n'+'time='+str(slp.time.values)+'\n'+'ABL'
        #plt.title(title)
        #plt.axis("off")
        ax0 = fig.add_subplot(gs[0,0])
        cs = ax0.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='bwr',vmax=290,vmin=250)
        cs0 = ax0.contour(x1,y1,znew0, levels=np.arange(250,290,5),colors='k')
        ax0.clabel(cs0, inline_spacing=0, fontsize=10, fmt="%.0f")
        ax0.annotate('a)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        ax0.set_title('ERA5')
        zipped_lists = zip(x, y,np.array(t850_section))
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear')
        ax1 = fig.add_subplot(gs[0,1])
        cs = ax1.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='bwr',vmax=290,vmin=250)
        cs1 = ax1.contour(x1,y1,znew0, levels=np.arange(250,290,5),colors='k')
        ax1.clabel(cs1, inline_spacing=0, fontsize=10, fmt="%.0f")
        ax1.annotate('b)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        ax1.set_title(title)
        zipped_lists = zip(x_a500, y_a500,np.array(t850_section_a500))
        sorted_pairs = sorted(zipped_lists)
        print(np.array(t850_section_a500))
        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear')
        ax2 = fig.add_subplot(gs[0,2])
        cs2 = ax2.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='bwr',vmax=290,vmin=250)
        cs3 = ax2.contour(x1,y1,znew0, levels=np.arange(250,290,5),colors='k')
        ax2.clabel(cs3, inline_spacing=0, fontsize=10, fmt="%.0f")
        ax2.annotate('c)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        ax2.set_title('A500')
        zipped_lists = zip(x_at, y_at,np.array(t850_section_at))
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear')
        ax3 = fig.add_subplot(gs[0,3])
        cs3 = ax3.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='bwr',vmax=290,vmin=250)
        cs4 = ax3.contour(x1,y1,znew0, levels=np.arange(250,290,5),colors='k')
        ax3.clabel(cs4, inline_spacing=0, fontsize=10, fmt="%.0f")
        ax3.annotate('d)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        ax3.set_title('AT')
        #fig.colorbar(cs2,ax=ax2,label= '[K]',shrink=.5,extend='both')
        ax_color = fig.add_subplot(gs[0, 4]) 
        cbar = fig.colorbar(cs2, cax=ax_color, shrink=.5,extend='both')
        cbar.set_label('[K]')
        #ax_color.axis("off")
        zipped_lists = zip(x_era5, y_era5,tpw_section_era5)
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        ax3 = fig.add_subplot(gs[1,0])
        cs = ax3.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],vmin=5,vmax=45)
        ax3.annotate('e)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        zipped_lists = zip(x, y,tpw_section)
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        ax4 = fig.add_subplot(gs[1,1])
        cs = ax4.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],vmin=5,vmax=45)
        ax4.annotate('f)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        zipped_lists = zip(x_a500, y_a500,tpw_section_a500)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        ax5 = fig.add_subplot(gs[1,2])
        cs5 = ax5.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],vmin=5,vmax=45)
        ax5.annotate('g)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        zipped_lists = zip(x_at, y_at,tpw_section_at)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        ax6 = fig.add_subplot(gs[1,3])
        cs6 = ax6.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],vmin=5,vmax=45)
        ax_color2 = fig.add_subplot(gs[1, 4]) 
        #ax_color2.axis("off")
        cbar = fig.colorbar(cs6, cax=ax_color2, shrink=.5,extend='both')
        cbar.set_label(r'[kg/m$^2$]')
        ax6.annotate('h)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        #extend='both'r'$[W/m^2]$'
        save_fn =  '%s_[%s].png' % ("t850_tpw_north_cyc2_20_updated2",str(int((4*k+((a+6)/6)))).zfill(3)) #479
        print('hi')
        plt.savefig(save_fn)
        zipped_lists = zip(x_era5, y_era5,np.array(u10_section_era5))
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear')
        zipped_lists2 = zip(x_era5, y_era5,np.array(slp_v_era5))
        sorted_pairs2 = sorted(zipped_lists2)
        tuples2 = zip(*sorted_pairs2)
        list4, list5,list6 = [ list(tuple) for tuple in  tuples2]
        znew1 = griddata((list4, list5), list6, (xv, yv), method='linear',fill_value=np.nan)
        fig = plt.figure(figsize=[21,14])
        gs = gridspec.GridSpec(2,5, wspace= 0.5,hspace=None,width_ratios=[1,1,1,1,0.05])
        title='10 m Wind Speed and 500 hPa Vertical Velocity'+'\n'+'time='+str(slp.time.values)+'\n'+'FV ABL'
        #plt.title(title)
        #plt.axis("off")
        ax0 = fig.add_subplot(gs[0,0])
        cs = ax0.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='YlOrRd',vmin=0,vmax=20)
        cs0 = ax0.contour(x1,y1,znew1, levels=np.arange(970,1015,5),colors='k')
        ax0.clabel(cs0, inline_spacing=0, fontsize=10, fmt="%.0f")
        ax0.annotate('a)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        ax0.set_title('ERA5')
        zipped_lists = zip(x, y,np.array(u10_section))
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear')
        zipped_lists2 = zip(x, y,np.array(slp_v))
        sorted_pairs2 = sorted(zipped_lists2)
        tuples2 = zip(*sorted_pairs2)
        list4, list5,list6 = [ list(tuple) for tuple in  tuples2]
        znew1 = griddata((list4, list5), list6, (xv, yv), method='linear',fill_value=np.nan)
        ax1 = fig.add_subplot(gs[0,1])
        cs = ax1.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='YlOrRd',vmin=0,vmax=20)
        cs1 = ax1.contour(x1,y1,znew1, levels=np.arange(970,1015,5),colors='k')
        ax1.clabel(cs1, inline_spacing=0, fontsize=10, fmt="%.0f")
        ax1.annotate('b)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        ax1.set_title(title)
        zipped_lists = zip(x_a500, y_a500,np.array(u10_section_a500))
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear')
        zipped_lists2 = zip(x_a500, y_a500,np.array(slp_v_a500))
        sorted_pairs2 = sorted(zipped_lists2)
        tuples2 = zip(*sorted_pairs2)
        list4, list5,list6 = [ list(tuple) for tuple in  tuples2]
        znew1 = griddata((list4, list5), list6, (xv, yv), method='linear',fill_value=np.nan)
        ax2 = fig.add_subplot(gs[0,2])
        cs2 = ax2.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='YlOrRd',vmin=0,vmax=20)
        cs3 = ax2.contour(x1,y1,znew1, levels=np.arange(970,1015,5),colors='k')
        ax2.clabel(cs3, inline_spacing=0, fontsize=10, fmt="%.0f")
        ax2.annotate('c)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        ax2.set_title('A500')
        zipped_lists = zip(x_at, y_at,np.array(u10_section_at))
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear')
        zipped_lists2 = zip(x_at, y_at,np.array(slp_v_at))
        sorted_pairs2 = sorted(zipped_lists2)
        tuples2 = zip(*sorted_pairs2)
        list4, list5,list6 = [ list(tuple) for tuple in  tuples2]
        znew1 = griddata((list4, list5), list6, (xv, yv), method='linear',fill_value=np.nan)
        ax3 = fig.add_subplot(gs[0,3])
        cs3 = ax3.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='YlOrRd',vmin=0,vmax=20)
        cs4 = ax3.contour(x1,y1,znew1, levels=np.arange(970,1015,5),colors='k')
        ax3.clabel(cs4, inline_spacing=0, fontsize=10, fmt="%.0f")
        ax3.annotate('d)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        ax3.set_title('AT')
        #fig.colorbar(cs2,ax=ax2,label= '[K]',shrink=.5,extend='both')
        ax_color = fig.add_subplot(gs[0, 4]) 
        cbar = fig.colorbar(cs3, cax=ax_color, shrink=.5,extend='both')
        cbar.set_label('[m/s]')
        #ax_color.axis("off")
        zipped_lists = zip(x_era5, y_era5,o_section_era5)
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        ax4 = fig.add_subplot(gs[1,0])
        cs = ax4.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='PiYG',vmin=-30,vmax=30)
        ax4.annotate('e)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        zipped_lists = zip(x, y,o_section)
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        ax5 = fig.add_subplot(gs[1,1])
        cs = ax5.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='PiYG',vmin=-30,vmax=30)
        ax5.annotate('f)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        zipped_lists = zip(x_a500, y_a500,o_section_a500)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        ax6 = fig.add_subplot(gs[1,2])
        cs6 = ax6.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='PiYG',vmin=-30,vmax=30)
        ax6.annotate('g)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        zipped_lists = zip(x_at, y_at,o_section_at)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        ax7 = fig.add_subplot(gs[1,3])
        cs7 = ax7.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='PiYG',vmin=-30,vmax=30)
        ax_color2 = fig.add_subplot(gs[1, 4]) 
        #ax_color2.axis("off")
        cbar = fig.colorbar(cs7, cax=ax_color2, shrink=.5,extend='both')
        cbar.set_label('[hPa/hour]')
        ax7.annotate('h)', xy=(0, 1), xycoords='axes fraction', xytext=(+0.5, -0.5), textcoords='offset fontsize',fontsize='medium', verticalalignment='top', fontfamily='serif')
        #extend='both'r'$[W/m^2]$'
        save_fn =  '%s_[%s].png' % ("u10_o_north_cyc2_20_updated2",str(int((4*k+((a+6)/6)))).zfill(3)) #479
        print('hi')
        plt.savefig(save_fn)
        '''zipped_lists = zip(x, y,slp_v)
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        fig = plt.figure()
        cs = plt.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='plasma',vmin=980,vmax=1010)
        cbar = plt.colorbar(cs,label='[hPa]')#extend='both'r'$[W/m^2]$'
        title='Cyclone 2 Sea Level Pressure'+'\n'+'time='+str(slp.time.values)
        plt.title(title)
        save_fn =  '%s_[%s].png' % ("slp_north_cyc2_25_updated",str(int((4*k+((a+6)/6)))).zfill(3)) #479
        plt.savefig(save_fn)'''
        
        '''zipped_lists = zip(x, y,np.array(u10_section))
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        zipped_lists2 = zip(x, y,np.array(slp_v))
        sorted_pairs2 = sorted(zipped_lists2)
        tuples2 = zip(*sorted_pairs2)
        list4, list5,list6 = [ list(tuple) for tuple in  tuples2]
        znew1 = griddata((list4, list5), list6, (xv, yv), method='linear',fill_value=np.nan)
        fig = plt.figure()
        cs1 = plt.contour(x1,y1,znew1, levels=np.arange(970,1015,5),colors='k')
        plt.clabel(cs1, inline_spacing=0, fontsize=10, fmt="%.0f")
        cs = plt.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='YlOrRd',vmin=0,vmax=20)
        cbar = plt.colorbar(cs,label='[m/s]',extend='max')#extend='both'r'$[W/m^2]$'
        title='Cyclone 2 10 meter Wind Speed'+'\n'+'time='+str(slp.time.values)
        plt.title(title)
        save_fn =  '%s_[%s].png' % ("u10_north_cyc2_25_updated",str(int((4*k+((a+6)/6)))).zfill(3)) #479
        plt.savefig(save_fn)
        zipped_lists = zip(x, y,np.array(o_section))
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        zipped_lists2 = zip(x, y,np.array(slp_v))
        sorted_pairs2 = sorted(zipped_lists2)
        tuples2 = zip(*sorted_pairs2)
        list4, list5,list6 = [ list(tuple) for tuple in  tuples2]
        znew1 = griddata((list4, list5), list6, (xv, yv), method='linear',fill_value=np.nan)
        fig = plt.figure()
        #cs1 = plt.contour(x1,y1,znew1, levels=np.arange(970,1015,5),colors='k')
        #plt.clabel(cs1, inline_spacing=0, fontsize=10, fmt="%.0f")
        cs = plt.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='PiYG',vmin=-30,vmax=30)
        cbar = plt.colorbar(cs,label='[hPa/hour]',extend='both')#extend='both'r'$[W/m^2]$'
        cbar.ax.invert_yaxis()
        title='Cyclone 2 500 hPa Vertical Velocity'+'\n'+'time='+str(slp.time.values)
        plt.title(title)
        save_fn =  '%s_[%s].png' % ("o_north_cyc2_25_updated",str(int((4*k+((a+6)/6)))).zfill(3)) #479
        plt.savefig(save_fn)
        plt.close()
        zipped_lists = zip(x, y,np.array(z_section))
        sorted_pairs = sorted(zipped_lists)

        tuples = zip(*sorted_pairs)
        list1, list2,list3 = [ list(tuple) for tuple in  tuples]
        x1 = np.arange(-1200, 1201, 50)
        y1 = np.arange(-1200, 1201, 50)
        xv, yv = np.meshgrid(x1, y1)
        znew0 = griddata((list1, list2), list3, (xv, yv), method='linear',fill_value=np.nan)
        zipped_lists2 = zip(x, y,np.array(z_section))
        sorted_pairs2 = sorted(zipped_lists2)
        tuples2 = zip(*sorted_pairs2)
        list4, list5,list6 = [ list(tuple) for tuple in  tuples2]
        znew1 = griddata((list4, list5), list6, (xv, yv), method='linear',fill_value=np.nan)
        fig = plt.figure()
        cs1 = plt.contour(x1,y1,znew1, levels=np.arange(500,600,4),colors='k')
        plt.clabel(cs1, inline_spacing=0, fontsize=10, fmt="%.0f")
        cs = plt.imshow(znew0, origin='lower',extent=[x1[0], x1[-1], y1[0], y1[-1]],cmap='plasma',vmin=500,vmax=600)
        cbar = plt.colorbar(cs,label='[dm]',extend='both')#extend='both'r'$[W/m^2]$'
        cbar.ax.invert_yaxis()
        title='Cyclone 2 500 hPa Geopotential Height'+'\n'+'time='+str(slp.time.values)
        plt.title(title)
        save_fn =  '%s_[%s].png' % ("z_north_cyc2_25_updated",str(int((4*k+(a+1)))).zfill(3)) #479
        plt.savefig(save_fn)
        plt.close()'''
        
