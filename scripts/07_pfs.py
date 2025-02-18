#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jason Chiappa 2025
# Examines the raw 1-hr stage iv files through each event and gathers statistics on the 1 mm/hr precipitation features (PFs)
# WARNING: Takes several hours to run!


# In[ ]:


import os, warnings
import numpy as np
import pandas as pd
import xarray as xr
import datetime
from datetime import datetime, timedelta
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

# Range of events to run [change if resuming]
# Can be 'start', 'end' or event id integer
start,end = 'start','end'

# Main data directory (containing ere_database_exattrs_[sep_dist].csv)
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/ERE_analysis/v2/data/'

# Accumulation map array directory
arr_dir = 'E:/Research/EREs/data/array_output/'

# directory for lat/lon grids
latlon_dir = 'E:/Research/EREs/data/stageiv_latlon/'

# st4 1h data directory
path = 'E:/Research/EREs/data/stage_iv/01h/'


# In[ ]:


##############################################################################################
#### Leave the stuff below ####
##############################################################################################

# input database file name
dbfname = 'ere_database_prelim_'+str(sep_dist)+'km.csv'

# output database file name
outputfname = 'ere_database_pfs_'+str(sep_dist)+'km.csv'


# In[ ]:


# degree longitude to distance (km)
def dlon_km(event_lat):
    coord1 = (event_lat,0)
    coord2 = (event_lat,1)
    return(geodesic(coord1,coord2).km)

# delete unnecessary idx files produced by xarray
def wipe_idx(path):
    files = [file for file in os.listdir(path) if file.endswith('.idx')]
    for file in files:
        os.remove(path+file)


# In[ ]:


wipe_idx(path)


# In[ ]:


# list of attribute names
attrs = ['event','time','pf_length','lat1','lon1','lat2','lon2','azimuth','pf_points','pf_average']

data = {}
for attr in attrs:
    data[attr] = []


# In[ ]:


if resume:
    # load existing and save parameters to list to be added to
    df0 = pd.read_csv(datadir+outputfname)
    # df0['time'] = pd.to_datetime(df0['time'], format="%Y%m%d%H")
    # save parameters for events before run period
    df0 = df0[df0.event < start]
    
    for attr in attrs:
        data[attr] += df0[attr].to_list()

# event database
df = pd.read_csv(datadir+dbfname)
df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")
df['start_time'] = pd.to_datetime(df['start_time'], format="%Y-%m-%d %H:%M:%S")
df['accum_time'] = pd.to_datetime(df['accum_time'], format="%Y-%m-%d %H:%M:%S")

if start=='start':
    start = df['event'].to_list()[0]
if end=='end':
    end = df['event'].to_list()[-1]

# dataframe with events to analyze
df = df[(df.event >= start) & (df.event <= end)]


# In[ ]:


# define lists from dataframe
events = df.event.to_list()
starttimes = df.start_time.to_list()
endtimes = df.accum_time.to_list()
event_lats = df.lat.to_list()
event_lons = df.lon.to_list()
latvers = df.latversion.to_list()


# In[ ]:


# groups all contiguous points exceeding a threshold surounding the point of max exceedance
# input precipitation map array, origin point lat/lon, threshold in mm, search radius in pixels (4km)
def group(accum_array,maxlat,maxlon,threshold=1,searchradius=1):
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


# In[ ]:


# search for nearby PF if not connected to point of max exceedance
# return new dummy maxlat and maxlon
def search_pf(accum_array,maxlat,maxlon,threshold=1):
    dlon = dlon_km(maxlat)
    dlat = 111
    
    # lat/lon box radius for 1 pixel search (4 km grid spacing)
    latrad = (4/dlat) + (2/dlat)
    lonrad = (4/dlon) + (2/dlon)

    thresh_exceed = accum_array - threshold

    latexceed = lat[thresh_exceed > 0]
    lonexceed = lon[thresh_exceed > 0]
    accumexceed = accum_array[thresh_exceed > 0]
    
    # if fail
    maxlat2 = None
    maxlon2 = None
    
    # attempt search within 2 pixels and repeat for higher radii if no exceedance found
    for multi in range(1,11):
        max_indexs1 = np.where((latexceed>maxlat-multi*latrad) & (latexceed<maxlat+multi*latrad) &
                     (lonexceed>maxlon-multi*lonrad) & (lonexceed<maxlon+multi*lonrad))[0]
        max_indexs2 = np.where((latexceed>maxlat-(multi+1)*latrad) & (latexceed<maxlat+(multi+1)*latrad) &
                     (lonexceed>maxlon-(multi+1)*lonrad) & (lonexceed<maxlon+(multi+1)*lonrad))[0]
        
        delete = [i for i in range(len(max_indexs2)) if max_indexs2[i] in max_indexs1]
        near_idxs = np.delete(max_indexs2,delete)

        if len(np.nonzero(accumexceed[near_idxs])[0]) > 0:
            # use maximum point in case of multiple
            newmaxidx = near_idxs[np.argmax(accumexceed[near_idxs])]
            maxlat2 = latexceed[newmaxidx]
            maxlon2 = lonexceed[newmaxidx]
            break
        else:
            continue

    return(maxlat2,maxlon2)


# In[ ]:


# Examine PFs for each event
for event in events:
    i = events.index(event)
    stime = starttimes[i]
    etime = endtimes[i]
    lat = np.load(latlon_dir+'stageiv_lats_'+str(int(latvers[i]))+'.npy')
    lon = np.load(latlon_dir+'stageiv_lons_'+str(int(latvers[i]))+'.npy')
    
    # get timestamps for each hour in the event's duration
    tstamps = np.arange(stime+pd.Timedelta(1,unit='H'),etime+pd.Timedelta(1,unit='H'),pd.Timedelta(1,unit='H'))
    times = [tstamp.item().strftime('%Y%m%d%H') for tstamp in tstamps]
    
    for time in times:
        ds = xr.open_dataset(path+'st4_conus.'+time+'.01h.grb2')
        if 'tp' in list(ds.variables):
            maxlat = event_lats[i]
            maxlon = event_lons[i]
    
            tp = ds.tp.data

            grouplats,grouplons,groupvals = group(tp,maxlat,maxlon)
            
            if len(grouplats) == 0:
                # search for nearest 1 mm/hr PF and redefine starting point
                maxlat,maxlon = search_pf(tp,maxlat,maxlon,threshold=1)
                
                if maxlat==None:
                    length,lat1,lon1,lat2,lon2,fwd_azimuth,pf_avg,pf_pts = -999,-999,-999,-999,-999,-999,-999,-999
                else:
                    grouplats,grouplons,groupvals = group(tp,maxlat,maxlon)
                
            if len(grouplats) > 0:
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
                    
                # get average value of pf and number of grid points in pf
                pf_avg = np.nanmean(groupvals)
                pf_pts = len(groupvals)
            
        else:
            length,lat1,lon1,lat2,lon2,fwd_azimuth,pf_avg,pf_pts = -999,-999,-999,-999,-999,-999,-999,-999
            
        data['event'].append(event)
        data['time'].append(time)
        data['pf_length'].append(length)
        data['lat1'].append(lat1)
        data['lon1'].append(lon1)
        data['lat2'].append(lat2)
        data['lon2'].append(lon2)
        data['azimuth'].append(fwd_azimuth)
        data['pf_points'].append(pf_pts)
        data['pf_average'].append(pf_avg)
    
    df_pfs = pd.DataFrame(data)
    df_pfs.to_csv(datadir+outputfname,index=False)

    wipe_idx(path)
    
    print(str(event)+'/'+str(events[-1]),end='\r')

