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

cyclone_centers=pd.DataFrame(columns=['datetime','lat','lon','slp','label'])
#files = [i for i in glob.glob('/glade/derecho/scratch/ajacobs/archive/release-cesm2.2.1_FHIST_f09_g17_L32_no_nudge_2013_2022/atm/hist/release-cesm2.2.1_FHIST_f09_g17_L32_no_nudge_2013_2022.cam.h2.2013*.nc')]
files = [i for i in glob.glob('2*nudge_a500_cyclone_centers.nc')]
files.sort()

label_number=1
for i,file in enumerate(files):
    test=xr.open_dataset(file)


#cyclone_centers=xr.open_dataset('era5_fv_august_2022_cyclone_centers.nc').to_dataframe()
    

    #for a in range(4):
    for a in range(len(test.datetime.values)):
#for a in range(3673):
        
        
        update_label=True        
        label = 'cyclone_'+str(label_number).zfill(3)
        if not cyclone_centers.empty:
            previous_cyclones=cyclone_centers.loc[cyclone_centers['datetime'] == (test.datetime.values[a]-pd.Timedelta(timedelta(hours=6)))]
            if not previous_cyclones.empty:
                for k in range(min(previous_cyclones.index),max(previous_cyclones.index)+1,1):
                    if findDistance(previous_cyclones.lat[k],previous_cyclones.lon[k],test.lat.values[a],test.lon.values[a]) <=600:
                        label = previous_cyclones.label[k]
                        update_label=False
                        break;
        cyclone_centers.loc[len(cyclone_centers)] =[test.datetime.values[a],test.lat.values[a],test.lon.values[a],test.slp.values[a],label]
        cyclone=cyclone_centers.to_xarray()
        cyclone.to_netcdf('all_nudge_a500_cyclone_centers.nc')
    
        if update_label:
            label_number=label_number+1
                
            
