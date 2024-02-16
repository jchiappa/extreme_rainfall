#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Jason Chiappa 2024
# Applies conditions for classification as tropcial cyclone, isolated, MCS/synoptic, and nocturnal
# Also adds diurnal cycle attributes and the following attributes based on each event's 1 mm/hr precipitation feature:
# Maximum major axis length, number of consicutive hours with length (km) > MCS/isolated threshold


# In[2]:


import os, warnings
import numpy as np
import pandas as pd
from sun import sun
from datetime import datetime, timedelta
from itertools import groupby

warnings.simplefilter('ignore')


# In[3]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# Main data directory
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/Publication/ERE_analysis/data/'

# Will run for entire dataset (regardless of resuming) due to short runtime

##############################################################################################

# Tropical cyclone center within tc_latlondist degrees (lat/lon) of point of max exceedance 
# and within tc_timerange hours of peak accumulation hour classified as TC
tc_latlondist = 3
tc_timerange = 24

# Isolated classified as having maximum PF length < iso_length (km)
iso_length = 200

# MCS/synoptic classified as having PF length >= mcs_length (km) for at least mcs_time consecutive hours
mcs_length = 200
mcs_time = 4

# Nocturnal classified as having peak accumulation hour at night and/or 
# having at least this fraction of nighttime hours to the total event duration
night_ratio = 0.5


# In[ ]:





# In[4]:


# open main ERE database
df = pd.read_csv(datadir+'ere_database_exceed_attrs.csv')

df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")
df['start_time'] = pd.to_datetime(df['start_time'], format="%Y-%m-%d %H:%M:%S")
df['accum_time'] = pd.to_datetime(df['accum_time'], format="%Y-%m-%d %H:%M:%S")

df


# In[5]:


events = df.event.to_list()
event_times = df.event_time.to_list()
start_times = df.start_time.to_list()
end_times = df.accum_time.to_list()
event_lats = df.max_lat.to_list()
event_lons = df.max_lon.to_list()


# In[ ]:





# In[6]:


# get number of hours between sunset and sunrise (location/date-dependent)
night_hours = []
# and flag whether peak accumulation hour is at night
night_flags = []

for i in range(len(events)):
    event_time = event_times[i]
    start_time = start_times[i]
    end_time = end_times[i]
    event_lat = event_lats[i]
    event_lon = event_lons[i]
    
    try:
        s = sun(lat=event_lat,long=event_lon)
        sset = s.sunset(datetime(event_time.year, event_time.month, event_time.day))
    except:
        s = sun(lat=event_lat,long=event_lon+360)
        sset = s.sunset(datetime(event_time.year, event_time.month, event_time.day))

    try:
        s = sun(lat=event_lat,long=event_lon)
        srise = s.sunrise(datetime(event_time.year, event_time.month, event_time.day))
    except:
        s = sun(lat=event_lat,long=event_lon+360)
        srise = s.sunrise(datetime(event_time.year, event_time.month, event_time.day))
        
    # check if event hour is at night
    h = event_time.hour
    
    if sset.hour>12:
        if ((h>sset.hour) | (h<srise.hour)):
            night_flags.append(1)
        else:
            night_flags.append(0)
    else:
        if ((h>sset.hour) & (h<srise.hour)):
            night_flags.append(1)
        else:
            night_flags.append(0)
    
    # count FULL hours between sunset and sunrise
    # nighttime hours are considered fully within the bounds of sunset-sunrise
    
    event_hours = pd.date_range(start_time,end_time-pd.Timedelta(1,unit='h'),freq='h').hour.to_list()
    
    if sset.hour>12:
        total_night_hours = sum((h>sset.hour) | (h<srise.hour) for h in event_hours)
    else:
        total_night_hours = sum((h>sset.hour) & (h<srise.hour) for h in event_hours)
        
    night_hours.append(total_night_hours)


# In[7]:


# add diurnal attributes
df2 = df.copy()
df2['night_flag'] = night_flags
df2['night_hours'] = night_hours


# In[8]:


# create nocturnal flags
noct_flags = []

for event in df2.event:
    event_df = df2[df2.event == event]
    
    if ((event_df.night_hours.to_list()[0]/event_df.event_duration.to_list()[0] >= night_ratio) | 
        (event_df.night_flag.to_list()[0] == 1)):
        noct_flags.append(1)
    else:
        noct_flags.append(0)


# In[ ]:





# In[9]:


# open tropical cyclone data
tc = pd.read_csv(datadir+'ibtracs.since1980.list.v04r00.csv', header=0, na_filter=False, low_memory=False)
tc = tc.drop(0, axis=0)

#Contert 'ISO_time' to datetime
tc['ISO_TIME'] = pd.to_datetime(tc['ISO_TIME'], format='%Y-%m-%d %H:%M:%S')

tc['LAT'] = pd.to_numeric(tc['LAT'])
tc['LON'] = pd.to_numeric(tc['LON'])

tc


# In[10]:


# get tc classifications
tc_flags = []

for i in range(len(df)):
    event = df.iloc[i]
    
    start = event.event_time - timedelta(hours=tc_timerange)
    end = event.event_time + timedelta(hours=tc_timerange)
    lat = event.max_lat
    lon = event.max_lon
    
    tc_hits = tc[(tc.ISO_TIME>=start) & (tc.ISO_TIME<=end) & (tc.LAT>=lat-tc_latlondist) & (tc.LAT<=lat+tc_latlondist) & 
             (tc.LON>=lon-tc_latlondist) & (tc.LON<=lon+tc_latlondist)]

    if len(tc_hits) > 0:
        tc_flags.append(1)
    else:
        tc_flags.append(0)


# In[ ]:





# In[11]:


# open PF database
df_pfs = pd.read_csv(datadir+'ere_database_pfs.csv')
df_pfs


# In[12]:


# get pf attributes
max_lengths = []
max_consecutive_durations = []
hour_counts = []
azimuths = []
# hours greater than iso length

for event in events:
    event_df = df_pfs[df_pfs.event == event]
    max_lengths.append(event_df.pf_length.max())
    
    lexceed_df = event_df[event_df.pf_length > mcs_length]
    
    hr_count = len(lexceed_df)
    hour_counts.append(hr_count)
    
    if hr_count == 1:
        azimuths.append(lexceed_df.azimuth.mean())
        max_consecutive_durations.append(1)

    elif hr_count > 1:
        azmax = lexceed_df.azimuth[lexceed_df.idxmax(axis=0)['pf_length']]
        azimuths.append(azmax)
        
        # find maximum number of consecutive hours exceeding threshold
        dlist = list(np.diff(lexceed_df.index))
        iterlist = [(k, sum(1 for i in g)) for k,g in groupby(dlist)]
        ilist2 = [b for b in iterlist if b[0]==1]
        if len(ilist2) > 0:
            ilist3 = [c[1] for c in ilist2]
            maxdur = max(ilist3)+1
            max_consecutive_durations.append(maxdur)
        else:
            max_consecutive_durations.append(1)
                        
    else:
        azimuths.append(-999)
        max_consecutive_durations.append(0)


# In[13]:


# add pf attributes
df2['max_pflength'] = max_lengths
# df2['length_hours'] = hour_counts
df2['pflength_duration'] = max_consecutive_durations
# df2['azimuth'] = azimuths


# In[ ]:





# In[14]:


# get isolated and mcs/synoptic classifications
iso_flags = []
mcs_flags = []

for event in events:
    event_df = df2[df2.event == event]
    
    if event_df.max_pflength.to_list()[0] < iso_length:
        iso_flags.append(1)
    else:
        iso_flags.append(0)
    
    if event_df.pflength_duration.to_list()[0] >= mcs_time:
        mcs_flags.append(1)
    else:
        mcs_flags.append(0)


# In[ ]:





# In[15]:


# add flags to dataframe
df_flagged = df2.copy()
df_flagged['tc_flag'] = tc_flags
df_flagged['iso_flag'] = iso_flags
df_flagged['mcs_flag'] = mcs_flags
df_flagged['noct_flag'] = noct_flags


# In[16]:


# edit iso/mcs flags so tc events do not conflict
iso_flags2 = []
mcs_flags2 = []

for event in events:
    event_df = df_flagged[df_flagged.event == event]
    
    if (event_df.iso_flag.to_list()[0]==1) & (event_df.tc_flag.to_list()[0]==0):
        iso_flags2.append(1)
    else:
        iso_flags2.append(0)
        
    if (event_df.mcs_flag.to_list()[0]==1) & (event_df.tc_flag.to_list()[0]==0):
        mcs_flags2.append(1)
    else:
        mcs_flags2.append(0)
        
df_flagged['iso_flag'] = iso_flags2
df_flagged['mcs_flag'] = mcs_flags2


# In[ ]:





# In[17]:


df_flagged


# In[ ]:





# In[18]:


df_flagged.to_csv(datadir+'ere_database_classflagged.csv',index=False)


# In[ ]:





# In[19]:


df_flagged.keys()


# In[ ]:




