#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jason Chiappa 2024
# Calculates the following attributes based on each event's exceedance map:
# Total exceedance volume, approximate major axis length of the exceedance swath and it's forward azimuth


# In[1]:


import os, warnings
import numpy as np
import pandas as pd
from geopy.distance import geodesic
from datetime import datetime
import pyproj
geodesicp = pyproj.Geod(ellps='WGS84')

warnings.simplefilter('ignore')


# In[6]:


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

# Accumulation map array directory
arr_dir = 'E:/Research/EREs/data/map_arrays/'

##############################################################################################

latlondir = datadir+'stageiv_latlon/'

# exceed array directory
path = arr_dir+'exceed_arrays/10yr/'


# In[3]:


# degree longitude to distance (km)
def dlon_km(event_lat):
    coord1 = (event_lat,0)
    coord2 = (event_lat,1)
    return(geodesic(coord1,coord2).km)


# In[4]:


# list of attribute names
attrs = ['ex_vol','ex_length','ex_azimuth']

data = {}
for attr in attrs:
    data[attr] = []

if resume:
    df = pd.read_csv(datadir+'ere_database_exceed_attrs.csv')
    df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")
    # save parameters for events before run period
    df = df[df.event_time < datetime(startyear,startmonth,1)]
    for attr in attrs:
        data[attr] += df[attr].to_list()

# raw dataframe
df0 = pd.read_csv(datadir+'ere_database_thresh.csv')
df0['event_time'] = pd.to_datetime(df0['event_time'], format="%Y-%m-%d %H:%M:%S")
# dataframe with events to analyze
df = df0[(df0.event_time >= datetime(startyear,startmonth,1)) & 
         (df0.event_time < datetime(endyear,endmonth,finalday) + timedelta(days=1))]
df


# In[5]:


events = df.event.to_list()
event_lats = df.max_lat.to_list()
event_lons = df.max_lon.to_list()
latvers = df.latversion.to_list()


# In[6]:


# groups all contiguous points exceeding a threshold surounding the point of max exceedance
# input precipitation map array, origin point lat/lon, threshold in mm, search radius in pixels (4km)
def group(accum_array,maxlat,maxlon,threshold=0,searchradius=1):
    dlon = dlon_km(maxlat)
    dlat = 111
    
    # lat/lon box radius for gathering adjacent/surrounding searchradius grid points only
    latrad = (searchradius*4/dlat) + (2/dlat)
    lonrad = (searchradius*4/dlon) + (2/dlon)

    thresh_exceed = accum_array - threshold

    latexceed = lat[thresh_exceed > 0]
    lonexceed = lon[thresh_exceed > 0]
    accumexceed = accum_array[thresh_exceed > 0]
    
    max_indexs = np.where((latexceed>maxlat-4/dlat/2) & (latexceed<maxlat+4/dlat/2) &
                         (lonexceed>maxlon-4/dlon/2) & (lonexceed<maxlon+4/dlon/2))[0]
    
    if len(max_indexs) == 0:
        return(np.array([]),np.array([]),np.array([]))
    
    else:
        max_index = max_indexs[0]

        grouplats = [latexceed[max_index]]
        grouplons = [lonexceed[max_index]]
        groupvals = [accumexceed[max_index]]

        latexceed = np.delete(latexceed,max_index)
        lonexceed = np.delete(lonexceed,max_index)
        accumexceed = np.delete(accumexceed,max_index)

        grouplen = [0,len(grouplats)]

        while grouplen[-1] > grouplen[-2]:
            for ptidx in range(grouplen[-2],grouplen[-1]):
                ptlat = grouplats[ptidx]
                ptlon = grouplons[ptidx]
                nearidxs = np.where((latexceed>ptlat-latrad) & (latexceed<ptlat+latrad) &
                                    (lonexceed>ptlon-lonrad) & (lonexceed<ptlon+lonrad))[0]
                if len(nearidxs) > 0:
                    for nidx in nearidxs:
                        grouplats.append(latexceed[nidx])
                        grouplons.append(lonexceed[nidx])
                        groupvals.append(accumexceed[nidx])
                latexceed = np.delete(latexceed,nearidxs)
                lonexceed = np.delete(lonexceed,nearidxs)
                accumexceed = np.delete(accumexceed,nearidxs)

            grouplen.append(len(grouplats))

        return(np.array(grouplats),np.array(grouplons),np.array(groupvals))


# In[7]:


def extreme_length(event):
    global lat
    global lon
    global tp
    
    i = events.index(event)
    maxlat = event_lats[i]
    maxlon = event_lons[i]
    lat = np.load(latlondir+'stageiv_lats_'+str(int(latvers[i]))+'.npy')
    lon = np.load(latlondir+'stageiv_lons_'+str(int(latvers[i]))+'.npy')

    tp = np.load(path+str(event).zfill(5)+'.npy')

    grouplats,grouplons,groupvals = group(tp,maxlat,maxlon)

    if len(grouplats) == 0:
        length,lat1,lon1,lat2,lon2,fwd_azimuth = -999,-999,-999,-999,-999,-999

    else:
        lat_argmin = np.nanargmin(grouplats)
        lat_argmax = np.nanargmax(grouplats)
        lon_argmin = np.nanargmin(grouplons)
        lon_argmax = np.nanargmax(grouplons)

        lat_minpt = (grouplats[lat_argmin],grouplons[lat_argmin])
        lat_maxpt = (grouplats[lat_argmax],grouplons[lat_argmax])
        lon_minpt = (grouplats[lon_argmin],grouplons[lon_argmin])
        lon_maxpt = (grouplats[lon_argmax],grouplons[lon_argmax])

        lat_dist = geodesic(lat_minpt,lat_maxpt).km
        lon_dist = geodesic(lon_minpt,lon_maxpt).km

        imax = [lat_dist,lon_dist].index(max([lat_dist,lon_dist]))
        length = [lat_dist,lon_dist][imax]
        lat1 = [lat_minpt[0],lon_minpt[0]][imax]
        lon1 = [lat_minpt[1],lon_minpt[1]][imax]
        lat2 = [lat_maxpt[0],lon_maxpt[0]][imax]
        lon2 = [lat_maxpt[1],lon_maxpt[1]][imax]

        fwd_azimuth,back_azimuth,distance = geodesicp.inv(lon1, lat1, lon2, lat2)

        if fwd_azimuth < 0:
            fwd_azimuth = fwd_azimuth+180
            
    if length==0.0:
        fwd_azimuth = -999
            
    return(length,fwd_azimuth)


# In[8]:


# perform calculations for each event
for event in events:
    exceed_array = np.load(path+str(event).zfill(5)+'.npy')
    
    # calculate exceedance volume
    data['ex_vol'].append(np.nansum(exceed_array))
    
    # calculate exceedance swath length and azimuth
    length,fwd_azimuth = extreme_length(event)
    data['ex_length'].append(length)
    data['ex_azimuth'].append(fwd_azimuth)
    
    print(str(event)+'/'+str(events[-1]),end='\r')


# In[9]:


df0


# In[10]:


for attr in attrs:
    df0[attr] = data[attr]

df0


# In[13]:


df0.to_csv(datadir+'ere_database_exceed_attrs.csv',index=False)


# In[ ]:




