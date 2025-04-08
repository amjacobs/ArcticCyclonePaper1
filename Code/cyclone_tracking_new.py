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
from cartopy.feature import NaturalEarthFeature
from netCDF4 import Dataset 
from scipy.stats import binned_statistic
import os
os.environ['ESMFMKFILE']='/glade/work/ajacobs/conda-envs/asp2023/lib/esmf.mk'
import xesmf as xe
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units
def findDistance(lat1,lon1,lat2,lon2):
  
    center = (lat1,lon1)
    point = (lat2,lon2)
    dist = distance.great_circle(center,point).km

    return dist
def process_slp_section(distance_array,slp,grad):
    distance_matrix = np.reshape(np.array(distance_array), (-1, 288))        
    slp_section=slp[np.where(distance_matrix<=1200)[0],np.where(distance_matrix<=1200)[1]].copy()
    if np.isnan(slp_section).any():
        grad = False 
        return [np.nan],[np.nan],slp,grad;
    grad_dist=distance_matrix[np.where(distance_matrix<=1200)[0],np.where(distance_matrix<=1200)[1]]
    slp[np.where(distance_matrix<=1200)[0],np.where(distance_matrix<=1200)[1]]=np.nan
    return slp_section,grad_dist,slp,grad;
#phis=xr.open_dataset('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_May_Sept2022_fv09/linear_full_nudge/atm/hist/fv9.lin.cam.h3.nc')
import warnings
warnings.filterwarnings("ignore")
cyclone_centers=pd.DataFrame(columns=['datetime','lat','lon','slp','label'])
files = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022.cam.h2.2013*.nc')]
#files = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h2.2013*.nc')]
files.sort()
phis_files = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_at_2013_2022.cam.h3.2013*.nc')]
#phis_files = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09/uvtnudge/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_nudge_2013_2022__fv09.cam.h3.2013*.nc')]
phis_files.sort()
label_number=1
for i,file in enumerate(files):
    phis=xr.open_dataset(phis_files[i])
    lats=np.where(phis.lat>=60)[0]
    test=xr.open_dataset(file)


#cyclone_centers=xr.open_dataset('era5_fv_august_2022_cyclone_centers.nc').to_dataframe()
    

    for a in range(4):
    #for a in [0,6,12,18]:
#for a in range(3673):
        grad=True
        geopotential_height=phis.Z3[a,-1,lats,:]

        geopotential_height=xr.where(geopotential_height>1500,0,1)
        slp=test.PSL[a,lats,:]
        i = 6
        while i > 0 and grad == True:
            update_label=True
            ind = np.where(slp==np.nanmin(slp))
            if len(ind[0])>1:
                ind=[np.array([ind[0][1]]),np.array([ind[1][1]])]
            if geopotential_height[ind[0],ind[1]]<1:
                slp[ind[0],ind[1]]=np.nan
            else:
                low=slp[ind[0],ind[1]].copy()
                slp_section = []
                grad_dist = []
                distance_array=[findDistance(low.lat,low.lon,slp[s][t].lat,slp[s][t].lon)
                        for s in range (len(slp))
                        for t in range(len(slp[0]))]
            
                slp_section,grad_dist,slp,grad=process_slp_section(distance_array,slp,grad)
    
                gradient=np.nanmean((np.array(slp_section)-low.values)*1000/np.array(grad_dist))/100
                if gradient >=7.5 and grad == True:
                    label = 'cyclone_'+str(label_number).zfill(3)
                    if not cyclone_centers.empty:
                        previous_cyclones=cyclone_centers.loc[cyclone_centers['datetime'] == (test.indexes['time'].values[a]-pd.Timedelta(timedelta(hours=6)))]
                        if not previous_cyclones.empty:
                            for k in range(min(previous_cyclones.index),max(previous_cyclones.index)+1,1):
                                if findDistance(previous_cyclones.lat[k],previous_cyclones.lon[k],low.lat.values[0],low.lon.values[0]) <=600:
                                    label = previous_cyclones.label[k]
                                    update_label=False
                                    break;
                    cyclone_centers.loc[len(cyclone_centers)] =[test.indexes['time'].values[a],low.lat.values[0],low.lon.values[0],low.values[0][0],label]
                    cyclone=cyclone_centers.to_xarray()
                    cyclone.to_netcdf('2013_nudge_a500_cyclone_centers.nc')
                    slp[ind[0],ind[1]] = np.nan
                    if update_label:
                        label_number=label_number+1
                else:
                    grad = False
            i =i-1
