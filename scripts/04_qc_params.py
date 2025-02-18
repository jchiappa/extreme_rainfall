#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jason Chiappa 2025
# Calculates parameters needed for quality control


# In[ ]:


import os, warnings
import numpy as np
import pandas as pd
import xarray as xr
from geopy.distance import geodesic
from datetime import datetime

warnings.simplefilter('ignore')


# In[ ]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# Separation distance between simultaneous events
sep_dist = 250
# sep_dist = 100

# IMPORTANT!!!
# True if resuming (or adding additional data to append to dataset), false if starting from beginning.
# Failing to change to True will overwrite all previous data and recalculate the parameters.
resume = False

# TIME RANGE for this run (YYYY-MM-DD HH:MM:SS) [change if resuming]
start_time = '2002-01-01 00:00:00'
end_time = '2024-12-31 23:00:00'

# Main data directory (containing ere_database_prelim_[sep_dist].csv)
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/ERE_analysis/v2/data/'

# Accumulation map array directory
arr_dir = 'E:/Research/EREs/data/array_output/'

# directory for lat/lon grids
latlon_dir = 'E:/Research/EREs/data/stageiv_latlon/'

# directory of US masks
usmask_dir = 'E:/Research/EREs/data/us_masks/'


# In[ ]:


##############################################################################################
#### Leave the stuff below ####
##############################################################################################

# 12-hr accumulation array subdirectory
accum_path = arr_dir+str(sep_dist)+'km/accum_12h/'

# preliminary database file name
dbfname = 'ere_database_prelim_'+str(sep_dist)+'km.csv'

# output database file name
outputfname = 'ere_database_qcparams_'+str(sep_dist)+'km.csv'

# list of all times (by hour) in selected time range
times = pd.Series(pd.date_range(start_time,end_time, freq='h')).to_list()

# lat/lon bounds for domain of study (W,E,S,N)
west,east,south,north = -104,-65,20,50


# In[ ]:


# degree longitude to distance (km)
def dlon_km(event_lat):
    coord1 = (event_lat,0)
    coord2 = (event_lat,1)
    return(geodesic(coord1,coord2).km)


# In[ ]:


if resume:
    # load existing and save parameters to list to be added to
    df0 = pd.read_csv(datadir+outputfname)
    df0['event_time'] = pd.to_datetime(df0['event_time'], format="%Y-%m-%d %H:%M:%S")
    # save parameters for events before run period
    df0 = df0[df0.event_time < times[0]]
    surround_avgs = df0.surround_avg.to_list()
    surround_mins = df0.surround_min.to_list()
else:
    # empty lists to store parameters
    surround_avgs = []
    surround_mins = []

# raw database
df = pd.read_csv(datadir+dbfname)
df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")

# dataframe with events to analyze
df1 = df[(df.event_time >= times[0]) & (df.event_time <= times[-1])]


# In[ ]:


# define lists from dataframe
events = df1.event.to_list()
event_lats = df1.lat.to_list()
event_lons = df1.lon.to_list()
latvers = df1.latversion.to_list()


# In[ ]:


# calculates parameters for quality control: 
# (1) average accumulation from st4 pixels surounding the point of max exceddance
# (2) minimum accumulation from st4 pixels surounding the point of max exceddance
# input event id number and search radius (number of pixels; keep as 1)
def qcparams(event,searchradius=1):
    i = events.index(event)
    lat = np.load(latlon_dir+'stageiv_lats_'+str(int(latvers[i]))+'.npy').flatten()
    lon = np.load(latlon_dir+'stageiv_lons_'+str(int(latvers[i]))+'.npy').flatten()
    
    # create lat/lon masks for domain
    latmask = (lat>=south) & (lat<=north)
    lonmask = (lon>=west) & (lon<=east)
    llmask = latmask & lonmask
    us_mask = np.load(usmask_dir+'US_mask_'+str(int(latvers[i]))+'.npy').flatten()
    mask = llmask.data & us_mask
    
    accum_array = np.load(accum_path+str(event).zfill(5)+'.npy').flatten()
    accum_array = np.where(mask, accum_array, np.nan)
    maxlat = event_lats[i]
    maxlon = event_lons[i]
    dlon = dlon_km(maxlat)
    dlat = 111
    max_index = np.where((lat>maxlat-4/dlat/2) & (lat<maxlat+4/dlat/2) &
                         (lon>maxlon-4/dlon/2) & (lon<maxlon+4/dlon/2))[0][0]
    # lat/lon box radius for gathering adjacent/surrounding searchradius grid points only
    latrad = (searchradius*4/dlat) + (2/dlat)
    lonrad = (searchradius*4/dlon) + (2/dlon)
    surround_idxs = np.where((lat>maxlat-latrad) & (lat<maxlat+latrad) &
                             (lon>maxlon-lonrad) & (lon<maxlon+lonrad))[0]
    surround_idxs = np.delete(surround_idxs,np.where(surround_idxs==max_index)[0][0])
    vals = accum_array[surround_idxs]
    event_surround_avg = np.nanmean(vals)
    event_surround_min = np.nanmin(vals)
    
    return(event_surround_avg,event_surround_min)


# In[ ]:


# perform calculation for each event
# long run time!
for event in events:
    event_surround_avg,event_surround_min = qcparams(event)
    surround_avgs.append(event_surround_avg)
    surround_mins.append(event_surround_min)
    print(str(event)+'/'+str(events[-1]),end='\r')


# In[ ]:


# add columns to original dataframe
df['surround_avg'] = surround_avgs
df['surround_min'] = surround_mins


# In[ ]:


# save to new csv
df.to_csv(datadir+outputfname,index=False)

