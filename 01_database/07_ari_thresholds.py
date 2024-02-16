#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jason Chiappa 2024
# Adds attributes for higher-end ARI thresholds at the points of max exceedance


# In[5]:


import os, warnings
import numpy as np
import pandas as pd
from datetime import datetime

warnings.simplefilter('ignore')


# In[2]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

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

# Main data directory
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/Publication/ERE_analysis/data/'

# NOAA Atlas 14 data directory
atlas_dir = 'E:/Research/EREs/data/noaa_atlas_14/'

##############################################################################################

latlondir = datadir+'stageiv_latlon/'


# In[17]:


# ARI thresholds to obtain
ARIs = ['25','50','100','500','1000']

thresh_data = {}
for ARI in ARIs:
    thresh_data['thresh'+ARI] = []
    
if resume:
    df = pd.read_csv(datadir+'ere_database_thresh.csv')
    df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")
    # save parameters for events before run period
    df = df[df.event_time < datetime(startyear,startmonth,1)]
    for ARI in ARIs:
        thresh_data['thresh'+ARI] += df['thresh'+ARI].to_list()

df0 = pd.read_csv(datadir+'ere_database_qc2.csv')
df0['event_time'] = pd.to_datetime(df0['event_time'], format="%Y-%m-%d %H:%M:%S")
df = df0[(df0.event_time >= datetime(startyear,startmonth,1)) & 
         (df0.event_time < datetime(endyear,endmonth,finalday) + timedelta(days=1))]

df


# In[18]:


events = df.event.to_list()
event_lats = df.max_lat.to_list()
event_lons = df.max_lon.to_list()
latvers = df.latversion.astype(int).to_list()


# In[19]:


# get thresholds at point of max exceedance for each event
for ARI in ARIs:
    print(ARI)
    thresh_list = []
    for i,event in enumerate(events):
        event_lat = event_lats[i]
        event_lon = event_lons[i]
        latver = latvers[i]
        
        if latver == 2:
            latver = 4
            
        lat = np.load(latlondir+'stageiv_lats_'+str(latver)+'.npy')
        lon = np.load(latlondir+'stageiv_lons_'+str(latver)+'.npy')
        
        atlas = np.load(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_'+str(ARI)+'yr_regrid_'+str(latver)+'.npy')
        
        ptidx = np.where((lat >= event_lat-1e-6) & (lat <= event_lat+1e-6) & (lon >= event_lon-1e-6) & (lon <= event_lon+1e-6))

        thresh_list.append(atlas[ptidx][0])
        print(str(i+1)+'/'+str(len(events)),end='\r')

    thresh_data['thresh'+ARI] += thresh_list
    df0['thresh'+ARI] = thresh_data['thresh'+ARI]
    print('\n')


# In[21]:


df0


# In[23]:


df0.to_csv(datadir+'ere_database_thresh.csv',index=False)


# In[ ]:




