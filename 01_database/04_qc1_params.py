#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Jason Chiappa 2024
# Calculates parameters needed for quality control


# In[1]:


import os, warnings
import numpy as np
import pandas as pd
import xarray as xr
from geopy.distance import geodesic
from datetime import datetime

warnings.simplefilter('ignore')


# In[3]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# IMPORTANT!!!
# True if resuming (or adding additional data to append to dataset), false if starting from beginning.
# Failing to change to true will overwrite all previous data and recalculate the parameters.
resume = False

# Enter start year and month for this run (change if resuming)
startyear = 2002
startmonth = 1

# Enter end year and month
endyear = 2023
endmonth = 12
# Last day of the month
finalday = 31

# Main data directory (containing ere_database_prelim.csv)
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/Publication/ERE_analysis/data/'

# Accumulation map array directory
arr_dir = 'E:/Research/EREs/data/map_arrays/'

##############################################################################################

latlondir = datadir+'stageiv_latlon/'

# 12-hr accumulation array subdirectory
accum_path = arr_dir+'accum_prelim_12h/'


# In[3]:


# degree longitude to distance (km)
def dlon_km(event_lat):
    coord1 = (event_lat,0)
    coord2 = (event_lat,1)
    return(geodesic(coord1,coord2).km)


# In[4]:


if resume:
    # load existing and save parameters to list to be added to
    df = pd.read_csv(datadir+'ere_database_prelim_qcparams.csv')
    df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")
    # save parameters for events before run period
    df0 = df[df.event_time < datetime(startyear,startmonth,1)]    
    surround_avgs = df0.surround_avg.to_list()
    surround_mins = df0.surround_min.to_list()
    
else:
    # empty lists to store parameters
    surround_avgs = []
    surround_mins = []

# raw dataframe
df1 = pd.read_csv(datadir+'ere_database_prelim.csv')
df1['event_time'] = pd.to_datetime(df1['event_time'], format="%Y-%m-%d %H:%M:%S")
# dataframe with events to analyze
df = df1[(df1.event_time >= datetime(startyear,startmonth,1)) & 
         (df1.event_time < datetime(endyear,endmonth,finalday) + timedelta(days=1))]
df


# In[5]:


# define lists from dataframe
events = df.event.to_list()
event_lats = df.max_lat.to_list()
event_lons = df.max_lon.to_list()
latvers = df.latversion.to_list()


# In[6]:


# calculates parameters for quality control: 
# (1) average accumulation from st4 pixels surounding the point of max exceddance
# (2) minimum accumulation from st4 pixels surounding the point of max exceddance
# input event id number and search radius (number of pixels; keep as 1)
def surround_avg(event,searchradius=1):
    i = events.index(event)
    lat = np.load(latlondir+'stageiv_lats_'+str(int(latvers[i]))+'.npy').flatten()
    lon = np.load(latlondir+'stageiv_lons_'+str(int(latvers[i]))+'.npy').flatten()
    accum_array =  np.load(accum_path+str(event).zfill(5)+'.npy').flatten()
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


# In[7]:


# perform calculation for each event
# long run time!
for event in events:
    event_surround_avg,event_surround_min = surround_avg(event)
    surround_avgs.append(event_surround_avg)
    surround_mins.append(event_surround_min)
    print(str(event)+'/'+str(events[-1]),end='\r')


# In[8]:


# add columns to raw dataframe
df1['surround_avg'] = surround_avgs
df1['surround_min'] = surround_mins


# In[9]:


# save to new csv
df1.to_csv(datadir+'ere_database_prelim_qcparams.csv',index=False)


# In[ ]:




