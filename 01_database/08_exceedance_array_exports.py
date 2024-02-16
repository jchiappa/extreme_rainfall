#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jason Chiappa 2024
# Exports map arrays of exceedance of 12-hr precip above the 10-yr ARI thresholds for each event


# In[5]:


import os, warnings
import numpy as np
import pandas as pd
from datetime import datetime

warnings.simplefilter('ignore')


# In[7]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

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

# Accumulation map array directory
arr_dir = 'E:/Research/EREs/data/map_arrays/'

##############################################################################################

latlondir = datadir+'stageiv_latlon/'
accum_path = arr_dir+'accum_12h/'
exceedpt_path = arr_dir+'exceedpts/10yr/'

# define export directories
exceed_dir = arr_dir+'exceed_arrays/10yr/'
os.makedirs(exceed_dir,exist_ok=True)


# In[9]:


# raw dataframe
df = pd.read_csv(datadir+'ere_database_thresh.csv')
df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")
# dataframe with events to analyze
df = df[(df.event_time >= datetime(startyear,startmonth,1)) & 
        (df.event_time < datetime(endyear,endmonth,finalday) + timedelta(days=1))]
df


# In[12]:


events = df.event.to_list()
latvers = df.latversion.astype(int).to_list()


# In[13]:


# generate exceedancy arrays for all 10-yr+ events
for event in events:
    i = events.index(event)
    
    lat = np.load(latlondir+'stageiv_lats_'+str(latvers[i])+'.npy')
    lon = np.load(latlondir+'stageiv_lons_'+str(latvers[i])+'.npy')
        
    # load dataset of 10-year recurrence thresholds (inches)
    atlas = np.load(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_10yr_regrid_'+str(latvers[i])+'.npy')

    # load 12-hr accumulation maps and exceedance points
    accum =  np.load(accum_path+str(event).zfill(5)+'.npy')
    exceed =  np.load(exceedpt_path+str(event).zfill(5)+'.npy')
    
    # mask out all points not included in this event (in case of simulataneous events)    
    mask = (lat==exceed[0][0]) & (lon==exceed[0][1])
    
    for expt in exceed[1:]:
        mask += (lat==expt[0]) & (lon==expt[1])
        
    exceed_array = mask*accum-atlas*mask
    
    np.save(exceed_dir+str(event).zfill(5)+'.npy',exceed_array)
    
    print(str(event)+'/'+str(events[-1]),end='\r')


# In[ ]:




