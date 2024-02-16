#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jason Chiappa 2024
# Combines events with SHARED EXCEEDANCE POINTS within 12 hours of each other.
# Events with higher maximum exceedance (retained) will have exceedance points counted outside of 250 km as well.
# Elimitated event accumulation maps are combined with retained map (still 12-hr accumulations) and exported to new array.
# Also exports numpy arrays with lists of exceedance point locations for each event.
# Since this completes the quality control, event ids are then redefined.
# WARNING: Takes several hours to run!


# In[1]:


import os, warnings
import shutil
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
from geopy.distance import geodesic

warnings.simplefilter('ignore')


# In[2]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# True if resuming (or adding additional data to append to dataset), false if starting from beginning.
# Failing to change to true will overwrite all previous data and recalculate the parameters.
resume = False

# Enter start year and month for this run (change if resuming)
# If resuming, recommend starting 1 month earlier in case events originally in different months can be combined
# Do not run on >1 month before end of preexisting dataset or event ids may get mixed up! (Must start over with resume = False)
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

# distance by which to separate simultaneous events (km) - i.e. radial distance from max accumulation (keep consistent)
sep_dist = 250

##############################################################################################

latlondir = datadir+'stageiv_latlon/'

# 12-hr accumulation array subdirectory
accum_path = arr_dir+'accum_prelim_12h/'
accum01_path = arr_dir+'accum_prelim_01h/'

# define export directories
accum_dir = arr_dir+'accum_12h/'
os.makedirs(accum_dir,exist_ok=True)
accum01_dir = arr_dir+'accum_01h/'
os.makedirs(accum01_dir,exist_ok=True)
exceedpt_dir = arr_dir+'exceedpts/10yr/'
os.makedirs(exceedpt_dir,exist_ok=True)

# temporary folders for exports before renaming to new event ids (in case of accidental replacement)
accum_dir_temp = arr_dir+'accum_12h_temp/'
os.makedirs(accum_dir_temp,exist_ok=True)
exceedpt_dir_temp = arr_dir+'exceedpts_temp/'
os.makedirs(exceedpt_dir_temp,exist_ok=True)


# In[3]:


if resume:
    # load existing and save parameters to list to be added to
    df = pd.read_csv(datadir+'ere_database_qc2_flagged.csv')
    df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")
    # save parameters for events before run period
    df0 = df[df.event_time < datetime(startyear,startmonth,1)]
    dup_flags = df0.dup_flag.to_list()
    max_sizes_new = df0.numpts.to_list()
    
else:
    dup_flags = []
    max_sizes_new = []
    # just in case folders have data, delete all files
    for f1,f2 in zip(os.listdir(accum_dir),os.listdir(exceedpt_dir)):
        os.remove(f1)
        os.remove(f2)

# raw dataframe
df1 = pd.read_csv(datadir+'ere_database_qc1.csv')
df1['event_time'] = pd.to_datetime(df1['event_time'], format="%Y-%m-%d %H:%M:%S")
# dataframe with events to analyze
df = df1[(df1.event_time >= datetime(startyear,startmonth,1)) & 
         (df1.event_time < datetime(endyear,endmonth,finalday) + timedelta(days=1))]
df


# In[4]:


# define lists from dataframe
events = df.event.to_list()
event_times = df.event_time.to_list()
event_lats = df.max_lat.to_list()
event_lons = df.max_lon.to_list()
max_exceeds = df.max_exceedance.to_list()
max_sizes = df.max_size.to_list()
etimesdti = pd.DatetimeIndex(event_times)
latvers = df.latversion.to_list()


# In[5]:


# get lat/lon points of exceedance for event (within sep_dist km)
def get_expoints(event):
    i = events.index(event)

    lat = np.load(latlondir+'stageiv_lats_'+str(int(latvers[i]))+'.npy')
    lon = np.load(latlondir+'stageiv_lons_'+str(int(latvers[i]))+'.npy')

    # load dataset of 10-year recurrence thresholds
    atlas = np.load(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_10yr_regrid_'+str(int(latvers[i]))+'.npy')

    maxlat = event_lats[i]
    maxlon = event_lons[i]

    # coordinate of max value
    maxcoord = (maxlat, maxlon)

    accum =  np.load(accum_path+str(event).zfill(5)+'.npy')

    # subtract threshold map
    exceed = accum-atlas

    # pull out points exceeding threshold
    extreme_mask = exceed>0

    exlats = lat[extreme_mask]
    exlons = lon[extreme_mask]

    eventpts = []
    # check each extreme point and save index if within sep_dist km
    for ptidx in range(len(exlats)):
        pt = (exlats[ptidx], exlons[ptidx])
        dist = geodesic(maxcoord, pt).km
        if dist < sep_dist:
            eventpts.append(ptidx)

    expts = [(exlats[eventpts][p],exlons[eventpts][p]) for p in range(len(exlats[eventpts]))]
    
    return(expts)


# In[6]:


def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
 
    if (a_set & b_set):
        common = list(a_set & b_set)
    else:
        common = []
        
    return(common)


# In[7]:


# list of events to flag as duplicates and remove
dup_events = []
# list of events to modify (combining with conflicting events)
combined_events = []
# list of new event exceedance point counts corresponding to combined_events
max_sizes_combined = []

# loop for each event
for event in events:
    
    i = events.index(event)

    maxexceed = max_exceeds[i]
    maxsize = max_sizes[i]

    accum =  np.load(accum_path+str(event).zfill(5)+'.npy')

    maxlat = event_lats[i]
    maxlon = event_lons[i]

    expts = get_expoints(event)

    # indexes of events to check (10 either side)
    irange = np.arange(i-10,i+10)[np.where((np.arange(i-10,i+10)>=0) & (np.arange(i-10,i+10)<len(events)) & 
                                           (np.arange(i-10,i+10)!=i))[0]]

    # list of event indexes that have matching points
    ematchidxs = []

    for i2 in irange:
        # time difference between events
        tdif = np.abs((etimesdti[i] - etimesdti[i2]) / np.timedelta64(1, 'h'))
        tcheck = tdif < 12

        # if within 12 hours check if matching exceedance points
        if tcheck:
            expts2 = get_expoints(events[i2])

            commonpts = common_member(expts, expts2)

            if len(commonpts) > 0:
                # save index of event
                ematchidxs.append(i2)

    # if any conflicting events, check if current event (i) has the maximum max exceedance value
    # if not, skip this event. If yes, combine accum arrays and save to new file, and update max_size parameter
    if len(ematchidxs) > 0:
        # check if exceedance/size of current event is greater than all other conflicting events
        exceedchecklist = []
        for matchidx in ematchidxs:
            maxexceed2 = max_exceeds[matchidx]
            exceedcheck = maxexceed > maxexceed2
            maxsize2 = max_sizes[matchidx]
            sizecheck = maxsize > maxsize2
            if exceedcheck | ((maxexceed >= maxexceed2) & sizecheck) :
                exceedchecklist.append(0) # if this event is highest, the sum of the list will be zero (no other event is greater)
            else:
                exceedchecklist.append(1)

        # case where it is the largest
        if sum(exceedchecklist) == 0:
            # save index to be retained
            combined_events.append(event)
            # combine accum arrays (keeping highest values only) and save all non-repeating exceedance points
            matchaccums = [np.nan_to_num(accum)]
            matchexceeds = [expts]
            for matchidx in ematchidxs:
                matchevent = events[matchidx]
                accummatch =  np.load(accum_path+str(matchevent).zfill(5)+'.npy')
                matchaccums.append(np.nan_to_num(accummatch))
                matchexceeds.append(get_expoints(matchevent))

            accum_combined = np.maximum.reduce(matchaccums)
            exceedpts_combined = np.array(list(set(sum(matchexceeds, [])))) # combines lists and removes repeating points
            max_sizes_combined.append(len(exceedpts_combined))

            np.save(accum_dir_temp+str(event).zfill(5)+'.npy',accum_combined)
            np.save(exceedpt_dir_temp+str(event).zfill(5)+'.npy',exceedpts_combined)

        # if event does not have a higher exceedance/size than all other conflicting events save index to discard
        else:
            dup_events.append(event)

    # else: keep event as is
    else:
        np.save(accum_dir_temp+str(event).zfill(5)+'.npy',np.nan_to_num(accum))
        np.save(exceedpt_dir_temp+str(event).zfill(5)+'.npy',np.array(expts))
        
    print(str(event)+'/'+str(events[-1]),end='\r')


# In[8]:


# create list of flags for events with shared exceedance (duplicates)
dup_flags += [int(e in dup_events) for e in events]  # 1 if duplicate event to be removed


# In[9]:


# new list of event sizes -- if event in combined_events, append new value from max_sizes_combined. else keep old value
for e in events:
    if e in combined_events:
        max_sizes_new.append(max_sizes_combined[combined_events.index(e)])
    else:
        max_sizes_new.append(max_sizes[events.index(e)])


# In[10]:


df1 = df1.drop('max_size',axis=1)
df1['dup_flag'] = dup_flags
df1['numpts'] = max_sizes_new
df1


# In[11]:


# save dataset with flags on duplicate events
df1.to_csv(datadir+'ere_database_qc2_flagged.csv',index=False)


# In[ ]:





# In[35]:


df_qc2 = df1[df1.dup_flag==0].drop(['dup_flag','surround_avg','surround_min'],axis=1)
df_qc2


# In[36]:


# renumber event ids
df_new = df_qc2.copy()

events_new = np.arange(1,len(df_new)+1)

df_new = df_new.reset_index(drop=True)
df_new['event'] = events_new
df_new


# In[37]:


# save dataset with event ids renumbered
df_new.to_csv(datadir+'ere_database_qc2.csv',index=False)


# In[32]:


# copy/rename necessary 1hr accumulation array files to new folder and rename
# move/rename new 12hr accumulation array and exceedance point files from temporary folders (run-time alert!)
events_qc2 = df_qc2[(df_qc2.event_time >= datetime(startyear,startmonth,1)) & 
                    (df_qc2.event_time < datetime(endyear,endmonth,finalday) + timedelta(days=1))].event.to_list()

events_new_files = df_new[(df_new.event_time >= datetime(startyear,startmonth,1)) & 
                          (df_new.event_time < datetime(endyear,endmonth,finalday) + timedelta(days=1))].event.to_list()

files = [str(f).zfill(5)+'.npy' for f in events_qc2]
files_new = [str(f).zfill(5)+'.npy' for f in events_new_files]

for file,file_new in zip(files,files_new):
    shutil.copyfile(accum01_path+file, accum01_dir+file_new)
    shutil.move(accum_dir_temp+file, accum_dir+file_new)
    shutil.move(exceedpt_dir_temp+file, exceedpt_dir+file_new)


# In[33]:


# delete temporary folders
shutil.rmtree(exceedpt_dir_temp)
shutil.rmtree(accum_dir_temp)


# In[ ]:




