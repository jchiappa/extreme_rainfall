#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jason Chiappa 2025
# Calculates the following attributes based on each event's exceedance footprint:
# Approximate major axis length of the exceedance swath associated with the point of max exceedance and it's forward azimuth


# In[ ]:


import os, warnings
import numpy as np
import pandas as pd
from geopy.distance import geodesic
import pyproj
geodesicp = pyproj.Geod(ellps='WGS84')

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


# In[ ]:


##############################################################################################
#### Leave the stuff below ####
##############################################################################################

# exceedance point list directory
path = arr_dir+str(sep_dist)+'km/exceedpts/10yr/'

# preliminary database file name
dbfname = 'ere_database_qcflagged_'+str(sep_dist)+'km.csv'

# output database file name
outputfname = 'ere_database_exattrs_'+str(sep_dist)+'km.csv'

# list of all times (by hour) in selected time range
times = pd.Series(pd.date_range(start_time,end_time, freq='h')).to_list()


# In[ ]:


# degree longitude to distance (km)
def dlon_km(event_lat):
    coord1 = (event_lat,0)
    coord2 = (event_lat,1)
    return(geodesic(coord1,coord2).km)


# In[ ]:


# list of attribute names
attrs = ['exceed_length','exceed_azimuth']

data = {}
for attr in attrs:
    data[attr] = []


# In[ ]:


if resume:
    # load existing and save parameters to list to be added to
    df0 = pd.read_csv(datadir+outputfname)
    df0['event_time'] = pd.to_datetime(df0['event_time'], format="%Y%m%d%H")
    # save parameters for events before run period
    df0 = df0[df0.event_time < times[0]]
    
    for attr in attrs:
        data[attr] += df0[attr].to_list()

# event database
df = pd.read_csv(datadir+dbfname)
df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")

# dataframe with events to analyze
df1 = df[(df.event_time >= times[0]) & (df.event_time <= times[-1])]


# In[ ]:


# define lists from dataframe
events = df1.event.to_list()
event_lats = df1.lat.to_list()
event_lons = df1.lon.to_list()


# In[ ]:


# groups all contiguous points exceeding a threshold surounding the point of max exceedance
# input event id, search radius in pixels (4km each)
def group(event,searchradius=1):
    i = events.index(event)

    maxlat = event_lats[i]
    maxlon = event_lons[i]
    
    dlon = dlon_km(maxlat)
    dlat = 111
    
    # lat/lon box radius for gathering adjacent/surrounding searchradius grid points only
    latrad = (searchradius*4/dlat) + (2/dlat)
    lonrad = (searchradius*4/dlon) + (2/dlon)

    expts = np.load(path+str(event).zfill(5)+'.npy')

    latexceed = expts[:,0]
    lonexceed = expts[:,1]
    
    max_indexs = np.where((latexceed>maxlat-4/dlat/2) & (latexceed<maxlat+4/dlat/2) &
                         (lonexceed>maxlon-4/dlon/2) & (lonexceed<maxlon+4/dlon/2))[0]
    

    max_index = max_indexs[0]

    grouplats = [latexceed[max_index]]
    grouplons = [lonexceed[max_index]]

    latexceed = np.delete(latexceed,max_index)
    lonexceed = np.delete(lonexceed,max_index)

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
            latexceed = np.delete(latexceed,nearidxs)
            lonexceed = np.delete(lonexceed,nearidxs)

        grouplen.append(len(grouplats))

        return(np.array(grouplats),np.array(grouplons))


# In[ ]:


def exceed_attrs(event):

    grouplats,grouplons = group(event)

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


# In[ ]:


# perform calculations for each event
for event in events:
    # calculate exceedance swath length and azimuth
    length,fwd_azimuth = exceed_attrs(event)
    data['exceed_length'].append(length)
    data['exceed_azimuth'].append(fwd_azimuth)
    
    print(str(event)+'/'+str(events[-1]),end='\r')


# In[ ]:


for attr in attrs:
    df[attr] = data[attr]


# In[ ]:


df.to_csv(datadir+outputfname,index=False)

