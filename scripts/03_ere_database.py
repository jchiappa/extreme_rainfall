#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jason Chiappa 2025
# This program creates a preliminary database of extreme precipitation events from the Stage IV dataset.
# EREs defined as a single or group of 12 hr accumulation grid points exceeding the 10-yr ARI threshold from NOAA Atlas 14.
# A csv file of event attributes, numpy arrays of maximum 12-hr precip, and accumulation map images are exported.
# Note: Takes around 2 hours per year of data

# 2025 version UPDATE:
# Uses an improved approach on gathering extreme events compared to the 2024 version (used for Chiappa et al. (2024) GRL paper)
# Fixes bug that causes certain simultaneous events/exceedance points to be excluded from the final database
# Improves on previous method by combining all events with shared exceedance points during this process instead of in post-processing
# This method helps prevent extreme 12-hr accumulations from being left out of the database
# Instead of just 1 specific event's 12-hr accumulation map, records the *maximum* 12-hr accumulation within a wider timeframe (when applicable)
# Adds parameters: earliest accumulation window start time, latest accumulation window end time
# Exports exceedance arrays and event map figures within program instead of during post-processing
# This version takes longer to run, but eliminates longer post-processing
# Precautions are in place for memory errors, but latest version makes them unlikely to occur

# Future research should implement the NOAA Atlas 15 dataset for thresholds when available


# In[ ]:


import os, warnings
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
from geopy.distance import geodesic
import geopandas as gpd
from shapely.geometry import Point, Polygon
from scipy import interpolate
from math import ceil, floor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import shapely.geometry as sgeom
from cartopy.geodesic import Geodesic
import gc
import psutil
from IPython.display import clear_output
import time as timemod

warnings.simplefilter('ignore')


# In[ ]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# Distance by which to separate simultaneous events (km) - i.e. radial distance from max accumulation
sep_dist = 250
# sep_dist = 100

# Original sep_dist = 250 km
# Recommended sep_dist = 100 km for future use
# This will create a larger database but will not merge simultaneous events with a relatively large separation.

# IMPORTANT!!!
# True if resuming (or adding additional data to append to dataset), false if starting from beginning.
# Failing to change to True will overwrite all previous data! Be extremely careful not to leave this as 'False'!
resume = False

# True if resuming after a memory error
resume_from_memory_error = False
# no need to change values below if True (will call correct values to resume from from memory error log)

# TIME RANGE for this run (YYYY-MM-DD HH:MM:SS) [see last entry in log file if debugging an error]
start_time = '2002-01-01 00:00:00'
end_time = '2024-12-31 23:00:00'

# Initialize lat/lon grid version id (the grid changes slightly through the dataset)
# latver = 1 if resume = False, else refer to last latver in previous run (see log file)
latver = 1

# Main data directory (for input and output of small files)
# Must include states_21basic folder (see README)
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/ERE_analysis/v2/data/'

# st4 1h data directory
path01 = 'E:/Research/EREs/data/stage_iv/01h/'

# st4 6h data directory
path06 = 'E:/Research/EREs/data/stage_iv/06h/'

# st4 24h data directory
path24 = 'E:/Research/EREs/data/stage_iv/24h/'

# NOAA Atlas 14 data directory
atlas_dir = 'E:/Research/EREs/data/noaa_atlas_14/'

# directory for lat/lon grids
latlon_dir = 'E:/Research/EREs/data/stageiv_latlon/'

# directory of US masks
usmask_dir = 'E:/Research/EREs/data/us_masks/'

# Map array export directory (will include folders for 12-hr accumulation maps and 1-hr maps at peak accumulation hours)
# Make sure to remove old version output or change the path first
arr_dir = 'E:/Research/EREs/data/array_output/'


# In[ ]:


##############################################################################################
#### Leave the stuff below ####
##############################################################################################

if resume_from_memory_error:
    start_time,latver = np.load(datadir+'mem_error_log.npy')
    start_time = str(start_time)
    latver = int(latver)

# directory for all array output
arr_dir += str(sep_dist)+'km/'

for fdir in [arr_dir,latlon_dir]:
    if (resume==False) & (os.path.exists(fdir)): # wipe all existing data!
        import shutil
        def delete_all_files_in_directory(directory):
            for filename in os.listdir(directory):
                file_path = os.path.join(directory, filename)
                try:
                    if os.path.isfile(file_path) or os.path.islink(file_path):
                        os.unlink(file_path)
                    elif os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                except Exception as e:
                    print(f'Failed to delete {file_path}. Reason: {e}')
        delete_all_files_in_directory(fdir)

# create directories if they do not exist
os.makedirs(arr_dir,exist_ok=True)
os.makedirs(latlon_dir,exist_ok=True)
os.makedirs(usmask_dir,exist_ok=True)
os.makedirs(atlas_dir,exist_ok=True)

# 12-hr accumulation map directory
accum_path = arr_dir+'accum_12h/'
os.makedirs(accum_path,exist_ok=True)

# exceedance map directory
exceed_path = arr_dir+'exceed_arrays/10yr/'
os.makedirs(exceed_path,exist_ok=True)

# exceedance point list directory
exceedpt_path = arr_dir+'exceedpts/10yr/'
os.makedirs(exceedpt_path,exist_ok=True)

# event map figure export directory
fig_path = arr_dir+'event_maps/'
os.makedirs(fig_path,exist_ok=True)

# log directory
log_dir = datadir+'logs/'
os.makedirs(log_dir,exist_ok=True)

# 1h st4 files
files = [f for f in os.listdir(path01) if f.endswith('.grb2')]

# 6h st4 files
files6 = [f for f in os.listdir(path06) if f.endswith('.grb2')]

# 24h st4 files
files24 = [f for f in os.listdir(path24) if f.endswith('.grb2')]

# create list of all times (by hour) to analyze
all_times = pd.Series(pd.date_range(start_time,end_time, freq='h'))

starttime = list(all_times)[0]
times = [t for t in all_times if ((t.month==starttime.month) & (t.year==starttime.year))]

# database file name
dbfname = 'ere_database_prelim_'+str(sep_dist)+'km.csv'

# lat/lon bounds for domain of study (W,E,S,N)
west,east,south,north = -104,-65,20,50


# In[ ]:


# open shapefile for creating CONUS masks
shapefile = gpd.read_file(datadir+'states_21basic/states.shp')
states = shapefile.geometry.to_list()
state_abbr = shapefile.STATE_ABBR.to_list()

# check if point is in one of the states then return state abreviation, else None
def in_us(lat, lon):
    p = Point(lon, lat)
    for i in range(1,50):
        if states[i].contains(p):
            return state_abbr[i]
    return None
    
# delete unnecessary idx files produced by xarray
def wipe_idx(path):
    files = [file for file in os.listdir(path) if file.endswith('.idx')]
    for file in files:
        os.remove(path+file)


# In[ ]:


# for retrieving initial lat/lon fields
ds = xr.open_dataset(path01+'st4_conus.'+times[0].strftime('%Y%m%d%H')+'.01h.grb2')
lat = ds.latitude.data
lon = ds.longitude.data - 360


# In[ ]:


# interpolate NOAA Atlas 14 data (12-hr precipitation frequency) onto stage iv grid
def atlas_regrid():
    global atlas_dir
    global latver
    global lat
    global lon
    
    # load noaa altas 14 dataset
    ds_atlas = xr.open_dataset(atlas_dir+'NOAA_Atlas_14_CONUS.nc')
    
    alat = ds_atlas.lat.data
    alon = ds_atlas.lon.data
    
    # get 12-hr data and convert inches to mm
    atlas_pf12 = ds_atlas.pf_012_hr*25.4
    ARIs = atlas_pf12.ari
    
    # export arrays for the following ARIs with pf_012_hr
    ARI_selections = [10,25,50,100,500,1000]

    for ARI in ARI_selections:
        atlas_data = atlas_pf12[:,:,ARIs==ARI]

        # replace nan with high value to prevent triggers at coastline
        atlas_nonan = np.nan_to_num(atlas_data.data,nan=1e3)

        # regrid/interpolate atlas to match stage-iv
        f = interpolate.interp2d(alon,alat,atlas_nonan,kind='linear')

        atlas_regrided = np.empty_like(lat)

        for latix in range(atlas_regrided.shape[0]):
            for lonix in range(atlas_regrided.shape[1]):
                atlas_regrided[latix,lonix] = f(lon[latix,lonix],lat[latix,lonix])

        np.save(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_'+str(ARI)+'yr_regrid_'+str(latver)+'.npy',atlas_regrided)
        
        # clear RAM
        gc.collect()


# In[ ]:


# create a location mask (within USA borders/coastlines) for the current lat/lon grid
def us_mask(lat,lon):
    global usmask_dir
    global latver
    
    os.makedirs(usmask_dir,exist_ok=True)

    binary = np.zeros(lat.shape)

    for lt in range(lat.shape[0]):
        for ln in range(lat.shape[1]):
            # check if in US
            if in_us(lat[lt,ln],lon[lt,ln]) != None:
                binary[lt,ln] = 1
                
    # location mask
    lmask = np.ma.make_mask(binary)
    
    np.save(usmask_dir+'US_mask_'+str(latver)+'.npy',lmask)


# In[ ]:


# export lat/lon grids if not already done for this lat/lon version... may take several minutes
if not os.path.exists(latlon_dir+'stageiv_lats_'+str(latver)+'.npy'):
    np.save(latlon_dir+'stageiv_lats_'+str(latver)+'.npy',lat)
    np.save(latlon_dir+'stageiv_lons_'+str(latver)+'.npy',lon)
    
    # regrid atlas and us mask data
    atlas_regrid()
    us_mask(lat,lon)


# In[ ]:


# create basic lat/lon masks for domain
latmask = (lat>=south) & (lat<=north)
lonmask = (lon>=west) & (lon<=east)
llmask = latmask & lonmask

# load regridded atlas arrays
atlas = {}
for ARI in [10,25,50,100,500,1000]:
    atlas[str(ARI)] = np.load(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_'+str(ARI)+'yr_regrid_'+str(latver)+'.npy')

# load the mask file
US_mask = np.load(usmask_dir+'US_mask_'+str(latver)+'.npy')

# combine lat/lon and US masks
mask = llmask.data & US_mask
mask_i = np.invert(mask)


# In[ ]:


# subroutine to retrieve 12hr precip from 1hr data files
# tp12 is the 12-hr accumulation array from the previous iteration
# tstr is the formated time string for the current iteration
def get_12hr_precip(tp12,tstr):
    global lat
    global lon
    # for first iteration, fill array with all 12 hours
    if np.nansum(tp12)==0:
        # fill array with data and sum 12 hours
        for idx,time12 in enumerate(times12):
            # time in string format
            tstr = time12.strftime('%Y%m%d%H')
            # 01h filename
            f = 'st4_conus.'+tstr+'.01h.grb2'
            # make sure data is not missing
            if os.path.exists(path01+f):
                # open data and get 1h precip array
                ds = xr.open_dataset(path01+f)
                lat = ds.latitude.data
                lon = ds.longitude.data - 360
                if 'tp' in list(ds.variables):
                    tp01 =  ds.tp.data
                else:
                    # fill with zeros
                    tp01 =  np.zeros((lat.shape[0], lat.shape[1]))
                    tp12[idx] = tp01
                ds.close()
                del(ds)
                # clear RAM
                gc.collect()
                
                tp12[idx] = tp01
            else:
                # fill with zeros
                tp01 =  np.zeros((lat.shape[0], lat.shape[1]))
                tp12[idx] = tp01

        tp = np.nansum(tp12,axis=0)
        
    else:
        # 01h filename
        f = 'st4_conus.'+tstr+'.01h.grb2'
        # make sure data is not missing
        if os.path.exists(path01+f):
            # open data and get 1h precip array
            ds = xr.open_dataset(path01+f)
            lat = ds.latitude.data
            lon = ds.longitude.data - 360
            if 'tp' in list(ds.variables):
                tp01 =  ds.tp.data
            else:
                # fill with zeros
                tp01 =  np.zeros((lat.shape[0], lat.shape[1]))
                # append new hour of data to previous array excluding the first hour
                tp12 = np.append(tp12[1:],np.array([tp01]),axis=0)
                tp = np.nansum(tp12,axis=0)
            ds.close()
            del(ds)
            # clear RAM
            gc.collect()
            
            # append new hour of data to previous array excluding the first hour
            tp12 = np.append(tp12[1:],np.array([tp01]),axis=0)
            tp = np.nansum(tp12,axis=0)
        else:
            # fill with zeros
            tp01 =  np.zeros((lat.shape[0], lat.shape[1]))
            # append new hour of data to previous array excluding the first hour
            tp12 = np.append(tp12[1:],np.array([tp01]),axis=0)
            tp = np.nansum(tp12,axis=0)
        
    return(tp12,tp,lat,lon)


# In[ ]:


# use to apply 6-hr data QC
def cross_check_06():
    global tp06sum
    # if same times as previous iteration, keep previous accumulation, else create new
    if not np.array_equal(times06,times06_old):
        # new precip array (from 6hr data)
        tp18 = np.empty((len(times06), lat.shape[0], lat.shape[1]))
        # fill array with data and sum
        for idx,time06 in enumerate(times06):
            # time in string format
            tstr = time06.strftime('%Y%m%d%H')
            # 06h filename
            f = 'st4_conus.'+tstr+'.06h.grb2'
            # make sure data is not missing
            if os.path.exists(path06+f):
                # open data and get 6h precip array
                ds = xr.open_dataset(path06+f)
                if 'tp' in list(ds.variables):
                    tp06 =  ds.tp.data
                    tp18[idx] = tp06
                else:
                    # fill with zeros
                    tp06 =  np.zeros((lat.shape[0], lat.shape[1]))
                    tp18[idx] = tp06
                ds.close()
                del(ds)
                # clear RAM
                gc.collect()
            else:
                # fill with zeros
                tp06 =  np.zeros((lat.shape[0], lat.shape[1]))
                tp18[idx] = tp06
            tp06sum = np.nansum(tp18,axis=0)
        
    tpdiff = accum_raw - tp06sum

    # subtract difference to replace all erroneously high values with 6-hr data values
    accum_qc06 = accum_raw - np.nan_to_num(tpdiff * (tpdiff > 0))
    
    return(accum_qc06)


# In[ ]:


# use to apply 24-hr data QC on 6-hr QC'd data
def cross_check_24():
    global tp24sum
    # if same times as previous iteration, keep previous accumulation, else create new
    if not np.array_equal(times24,times24_old):
        # new precip array (from 6hr data)
        tp48 = np.empty((len(times24), lat.shape[0], lat.shape[1]))
        # fill array with data and sum
        for idx,time24 in enumerate(times24):
            # time in string format
            tstr = time24.strftime('%Y%m%d')+'12'
            # 24h filename
            f = 'st4_conus.'+tstr+'.24h.grb2'
            # make sure data is not missing
            if os.path.exists(path24+f):
                # open data and get 6h precip array
                ds = xr.open_dataset(path24+f)
                if 'tp' in list(ds.variables):
                    tp24 =  ds.tp.data
                    tp48[idx] = tp24
                else:
                    # fill with zeros
                    tp24 =  np.zeros((lat.shape[0], lat.shape[1]))
                    tp48[idx] = tp24
                ds.close()
                del(ds)
                # clear RAM
                gc.collect()
            else:
                # fill with zeros
                tp24 =  np.zeros((lat.shape[0], lat.shape[1]))
                tp48[idx] = tp24
            tp24sum = np.nansum(tp48,axis=0)
        
    tpdiff = accum_qc06 - tp24sum

    # subtract difference to replace all erroneously high values with 24-hr data values
    accum = accum_qc06 - np.nan_to_num(tpdiff * (tpdiff > 0))
    
    return(accum)


# In[ ]:


def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
 
    if (a_set & b_set):
        common = list(a_set & b_set)
    else:
        common = []
        
    return(common)


# In[ ]:


# check sorted data for any events with shared exceedance points within accumulation window
# this happens when hourly events at the same location are not continuous (i.e. events occur more than once within accumulation period)
# merge those events and export new sorted data
def merge_repeats():
    # merged events database
    data = {}
    for col in cols[1:]:
        data[col] = []
    # array data
    adata = {}
    for arr in arrs:
        adata[arr] = []
        
    # list of indexes combined
    ematchidxall = []
    
    for i in range(len(data_sorted['lat'])):
        # if already been combined, skip
        if i in ematchidxall:
            continue
    
        # initialize list of matching event indexes
        ematchidxs = [i]
    
        # initialize list of exceedance points for event
        expts = arrays_sorted['exceed_points'][i].copy()
        
        # indexes of events to check for shared accumulation period
        check_next = 10
        # can increase this value if dealing with a very large number of events in a short timespan
        
        for i2 in [a for a in range(i+1,i+check_next+1) if ((a < len(data_sorted['lat'])) & (a not in ematchidxall))]:
            # check for shared exceedance points or max exceedance points within sep_dist
            if (data_sorted['period_start'][i2] < data_sorted['period_end'][i]) & (data_sorted['period_end'][i2] > data_sorted['period_start'][i]):
                commonpts = common_member(expts, arrays_sorted['exceed_points'][i2])

                # if events occur within 12 hours of each other
                event_time_condition = abs((data_sorted['event_time'][i] - data_sorted['event_time'][i2]).total_seconds()/3600) < 12

                # if max exceedance points are within sep_dist
                pt1 = (data_sorted['lat'][i], data_sorted['lon'][i])
                pt2 = (data_sorted['lat'][i2], data_sorted['lon'][i2])
                sep_dist_condition = geodesic(pt1, pt2).km < sep_dist
                
                if (len(commonpts) > 0) | (event_time_condition & sep_dist_condition):
                    # save index of event
                    ematchidxs.append(i2)
                    # add new exceedance points
                    expts += list(set(arrays_sorted['exceed_points'][i2]) - set(expts))
                    
        # if any matches, combine and save to new events database
        if len(ematchidxs) > 1:
            ematchidxall+=ematchidxs
            
            # combine accumulation arrays to save only max values at each gridpoint
            matchaccums = [np.nan_to_num(arrays_sorted['accum_arrays'][i]) for i in ematchidxs]
            accum_combined = np.maximum.reduce(matchaccums)
            
            # accumulation period start and end times
            period_start = np.nanmin([data_sorted['period_start'][i] for i in ematchidxs])
            period_end = np.nanmax([data_sorted['period_end'][i] for i in ematchidxs])
            
            # find idx with max exceedance
            matchexceeds = [data_sorted['exceedance'][i] for i in ematchidxs]
            emaxidx = ematchidxs[matchexceeds.index(max(matchexceeds))]
            
            # create exceedance arrays
            # mask out all points not included in this event (in case of simulataneous events)    
            exmask = (lat==expts[0][0]) & (lon==expts[0][1])
            for expt in expts[1:]:
                exmask += (lat==expt[0]) & (lon==expt[1])
            atlas_e = np.load(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_10yr_regrid_'+str(int(data_sorted['latversion'][emaxidx]))+'.npy')
            exceed_array = exmask*accum_combined-atlas_e*exmask
            
            # calculate exceedance volume
            exvolume = np.nansum(exceed_array)
            
            # record data for event
            data['lat'].append(data_sorted['lat'][emaxidx])
            data['lon'].append(data_sorted['lon'][emaxidx])
            data['state'].append(data_sorted['state'][emaxidx])
            data['event_time'].append(data_sorted['event_time'][emaxidx])
            data['start_time'].append(data_sorted['start_time'][emaxidx])
            data['accum_time'].append(data_sorted['accum_time'][emaxidx])
            data['duration'].append(data_sorted['duration'][emaxidx])
            data['period_start'].append(period_start)
            data['period_end'].append(period_end)
            data['exceedance'].append(data_sorted['exceedance'][emaxidx])
            data['accumulation'].append(data_sorted['accumulation'][emaxidx])
            data['exceed_pts'].append(len(expts))
            data['exceed_vol'].append(exvolume)
            data['max_1hr'].append(data_sorted['max_1hr'][emaxidx])
            data['ari_10'].append(data_sorted['ari_10'][emaxidx])
            data['ari_25'].append(data_sorted['ari_25'][emaxidx])
            data['ari_50'].append(data_sorted['ari_50'][emaxidx])
            data['ari_100'].append(data_sorted['ari_100'][emaxidx])
            data['ari_500'].append(data_sorted['ari_500'][emaxidx])
            data['ari_1000'].append(data_sorted['ari_1000'][emaxidx])
            data['latversion'].append(data_sorted['latversion'][emaxidx])
            
            # append arrays to lists
            adata['accum_arrays'].append(accum_combined)
            adata['exceed_arrays'].append(exceed_array)
            adata['exceed_points'].append(expts)

    # delete events that were combined and append merged event(s)
    for col in cols[1:]:
        for i in sorted(ematchidxall, reverse=True):
            del data_sorted[col][i]
        data_sorted[col] += data[col]
    
    for arr in arrs:
        for i in sorted(ematchidxall, reverse=True):
            del arrays_sorted[arr][i]
        arrays_sorted[arr] += adata[arr]

    # re-sort lists
    # sort each list by event time
    sortedIndex = sorted(range(len(data_sorted['event_time'])), key=lambda k: data_sorted['event_time'][k])
    data_sorted2 = {}
    for col in cols[1:]:
        data_sorted2[col] = [data_sorted[col][i] for i in sortedIndex]
    
    # for array lists
    arrays_sorted2 = {}
    for arr in arrs:
        arrays_sorted2[arr] = [arrays_sorted[arr][i] for i in sortedIndex]

    return(data_sorted2,arrays_sorted2)


# In[ ]:


# generate and export accumulation maps for each month
def generate_figures():
    # create figure
    fig = plt.figure(figsize=(8,5),dpi=300)
    
    events = data_sorted['event']
    for i in range(len(data_sorted['event'])):
        event = data_sorted['event'][i]
        print('Generating figures: '+str(event)+'/'+str(events[-1]),end='\r')
        accum = arrays_sorted['accum_arrays'][i]
        expts = arrays_sorted['exceed_points'][i]
        thresh = data_sorted['ari_10'][i]
        maxlat = data_sorted['lat'][i]
        maxlon = data_sorted['lon'][i]
        period_start = data_sorted['period_start'][i]
        period_end = data_sorted['period_end'][i]
        
        exlats = [l[0] for l in expts]
        exlons = [l[1] for l in expts]
        
        # create geo-axis
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
        
        # plot the data
        pcm = plt.pcolormesh(lon,lat,accum, cmap='cubehelix_r', transform=ccrs.PlateCarree(), 
                             vmin=0, vmax=thresh)
        
        plt.scatter(exlons,exlats,marker='o',facecolors='none',edgecolors='r',s=1.5,linewidths=0.4,
                    transform=ccrs.PlateCarree(),zorder=3)
        
        plt.title(str(event)+'   '+str(period_start)+' – '+str(period_end),loc='left')
        
        cbar = plt.colorbar(pcm,orientation='vertical',extend='max',pad=0.02)
        cbar.set_label(label='12-hr QPE (mm)')
        
        # add geography
        ax.coastlines(lw=0.5)
        ax.add_feature(cfeature.BORDERS,lw=0.5)
        ax.add_feature(cfeature.STATES,lw=0.5)
        
        extent = [min([(maxlon-4.5),(np.nanmin(exlons)-0.5)]),max([(maxlon+4.5),(np.nanmax(exlons)+0.5)]),
                  min([(maxlat-4.5),(np.nanmin(exlats)-0.5)]),max([(maxlat+4.5),(np.nanmax(exlats)+0.5)])]
        
        ax.set_extent(extent)
        
        gd = Geodesic()
        cp250 = gd.circle(lon=maxlon, lat=maxlat, radius=sep_dist*1000)
        geom = [sgeom.Polygon(cp250)]
        
        ax.add_geometries(geom, crs=ccrs.PlateCarree(), facecolor='none', edgecolor='k', alpha=0.5, linewidth = 0.5)
        
        ax.scatter([maxlon],[maxlat],marker='o',facecolors='none',edgecolors='white',s=1.5,linewidths=0.6,
                    transform=ccrs.PlateCarree(),zorder=3)
            
        plt.savefig(fig_path+str(event).zfill(5)+'.jpg', dpi=300, bbox_inches='tight')

        # clear figure contents so figure can be reused (saves RAM)
        plt.clf()
        
    plt.close()
    gc.collect()


# In[ ]:


# List of column names to record data
cols = ['event', 'lat', 'lon', 'state', 'event_time', 'start_time', 'accum_time', 'duration', 'period_start', 'period_end', 
        'exceedance', 'accumulation', 'exceed_pts', 'exceed_vol', 'max_1hr', 'ari_10', 'ari_25', 'ari_50', 'ari_100',
        'ari_500', 'ari_1000', 'latversion']
        
# PARAMETER DESCRIPTIONS
# 'event': event ID
# 'lat', 'lon', 'state': # lats/lons/state of point of max exceedance
# 'event_time': time with max hourly total (at point of max exceedance)
# 'start_time': start time of hourly precip >=1 mm at point of max exceedance
# 'accum_time': end time of 12-hr accumulation period at point of max exceedance
# 'duration': difference between start and accum time
# 'period_start': full accumulation period start time (first merged event accum_time -12 hrs)
# 'period_end': full accumulation period end times (last merged event accum_time)
# 'exceedance': max exceedance above the 10-yr ARI threshold
# 'accumulation': 12-hr accumulation value at point of max exceedance
# 'exceed_pts': total # of gridpoints exceeding 10-yr ARI threshold
# 'exceed_vol': sum of 10-yr ARI exceedance at all gridpoints in event
# 'max_1hr': maximum 1-hr rainfall at point of max exceedance
# 'ari_10', 'ari_25', 'ari_50', 'ari_100', 'ari_500', 'ari_1000': ARI thresholds at point of max exceedance for respective return periods
# 'latversion': lat/lon grid version ID (changes a few times through the stage iv dataset)


# In[ ]:


# CAREFUL: WILL REPACE PREVIOUS DATA. INITIALIZATION ONLY! (make sure resume = True if resuming)
# Creates new csv file with empty columns to store data
if not resume:
    # create dictionary of empty lists to store data
    data = {}
    for col in cols:
        data[col] = []

    # create dataframe
    df = pd.DataFrame(data)

    # SAVE TO NEW FILE
    df.to_csv(datadir+dbfname,index=False)
    
    # arrays to save for each event
    # merged/max 12-hr accumulation arrays (2d maps)
    accum_arrays = []
    # 2d maps of exceedance above 10-yr ARI thresholds
    exceed_arrays = []
    # lists of exceedance points
    exceed_points = []


# In[ ]:


##############################################################################################
#### Execute ####
##############################################################################################
code_start = timemod.time()
code_start_str = timemod.strftime("%Y%m%d%H%M%S", timemod.localtime())

current_time = timemod.strftime("%Y-%m-%d %H:%M:%S", timemod.localtime())

# create log file
with open(log_dir+'log_'+str(sep_dist)+'km_'+code_start_str+'.txt', 'w') as f:
    f.write('Start: '+current_time+'\n')

wipe_idx(path01)
wipe_idx(path06)
wipe_idx(path24)

# execute algorithm one month at a time until all times have been analyzed
while len(all_times)>0:
    elapsed_time = timemod.time() - code_start
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    time_elapsed = f"Elapsed time: {int(hours):02}:{int(minutes):02}:{int(seconds):02}"
    print(time_elapsed)
    
    # dataframe up to previous month
    df0 = pd.read_csv(datadir+dbfname)
    
    # Convert strings to datetime
    if len(df0)>0:
        df0['event_time'] = pd.to_datetime(df0['event_time'], format="%Y-%m-%d %H:%M:%S")
        df0['accum_time'] = pd.to_datetime(df0['accum_time'], format="%Y-%m-%d %H:%M:%S")
        df0['period_start'] = pd.to_datetime(df0['period_start'], format="%Y-%m-%d %H:%M:%S")
        df0['period_end'] = pd.to_datetime(df0['period_end'], format="%Y-%m-%d %H:%M:%S")

        # cut off previous database at start of period
        df0 = df0[df0['period_end']<all_times[0]]

        # Initialize event ID number
        event_id = df0['event'].to_list()[-1]+1

    else:
        event_id = 1
    
    # define lists for hourly preliminary events
    expts_h,maxpt_h,event_time_h,event_max_h,event_exceedance_h,maxrate_h,event_size_h,\
    event_state_h,start_time_h,accum_time_h,accum_array_h,latver_h = [],[],[],[],[],[],[],[],[],[],[],[]
    
    # define empty lists to store data (see above cell for list of columns and descriptions)
    data = {}
    for col in cols:
        data[col] = []
        
    # arrays to save for each event
    # merged/max 12-hr accumulation arrays (2d maps)
    accum_arrays = []
    # 2d maps of exceedance above 10-yr ARI thresholds
    exceed_arrays = []
    # lists of exceedance points
    exceed_points = []
    
    # initialize cross-check file times
    times06_old = 0
    times24_old = 0
    
    # 12hr precip array
    tp12 = np.empty((12, lat.shape[0], lat.shape[1]))
    
    # obtain list of timestamps for month to iterate through
    mstarttime = list(all_times)[0]
    times = [t for t in all_times if ((t.month==mstarttime.month) & (t.year==mstarttime.year))]
    
    # loop for each hour in the month
    for index,time in enumerate(times):
        # show progress percentage for the month
        print(str(time)+' '+"{:.2f}".format((index+1)/len(times)*100)+'% ',end='\r')
        
        # time in string format
        tstr = time.strftime('%Y%m%d%H')
        
        # list of times in accumulation period
        times12 = pd.Series(pd.date_range(time-pd.Timedelta(11,'h'),time, freq='h'))
        
        # list of times for 6hr data files -- up to 3 (note that file names are the END times of 6h periods)
        times06 = pd.Series(pd.date_range((time-pd.Timedelta(6,'h')).floor('6h'),
                                          (time+pd.Timedelta(5,'h')).floor('6h'), freq='6h'))
        # same for 24hr data files -- up to 2 (note that 24hr files run from 12z to 12z each day)
        times24 = pd.Series(pd.date_range(time.floor('24h'),(time+pd.Timedelta(11,'h')).floor('24h')))
    
        # previous latitude array
        lat_prev = lat.copy()
    
        # get data
        tp12,accum_raw,lat,lon = get_12hr_precip(tp12,tstr)
    
        # save new lat/lon arrays if lat array is different from previous
        if np.array_equal(lat,lat_prev) == False:
            latver+=1
            if not os.path.exists(latlon_dir+'stageiv_lats_'+str(latver)+'.npy'):
                print('\nGenerating new lat/lon version grids')
                np.save(latlon_dir+'stageiv_lats_'+str(latver)+'.npy',lat)
                np.save(latlon_dir+'stageiv_lons_'+str(latver)+'.npy',lon)
                # redo usmask
                us_mask(lat,lon)
                # regrid atlas data
                atlas_regrid()
            
            latmask = (lat>=south) & (lat<=north)
            lonmask = (lon>=west) & (lon<=east)
            llmask = latmask & lonmask
            
            US_mask = np.load(usmask_dir+'US_mask_'+str(latver)+'.npy')
            mask = llmask.data & US_mask
            mask_i = np.invert(mask)
            
            # create dictionary of atlas arrays
            atlas = {}
            for ARI in [10,25,50,100,500,1000]:
                atlas[str(ARI)] = np.load(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_'+str(ARI)+'yr_regrid_'+str(latver)+'.npy')
            
        # perform cross-checks on data to take advantage of manual QC on 6 and 24 hr data
        accum_qc06 = cross_check_06()
        accum = cross_check_24()        
        
        # save times for next iteration to avoid unnecessary repeated opening of files
        times06_old = times06.copy()
        times24_old = times24.copy()
    
        # apply mask
        accum[mask_i] = np.nan
    
        # subtract threshold map
        exceed = accum-atlas['10']
    
        # if any points exceeding threshold
        if np.nanmax(exceed) > 0:
            # pull out points exceeding threshold
            extreme_mask = exceed>0
            exlats = lat[extreme_mask]
            exlons = lon[extreme_mask]
            exvals = accum[extreme_mask]
            exexceed = exceed[extreme_mask]
    
            # check for simultaneous events within sep_dist km
            # list index of points in separate simultaneous events (all points for first iteration)
            separate = list(range(len(exvals)))
    
            # save original expoints
            exlats_orig = exlats.copy()
            exlons_orig = exlons.copy()
            expts_orig = [(exlats[i], exlons[i]) for i in range(len(exlats))]
    
            # while points not sampled yet, count additional events
            while len(separate) > 0:
                # index points that have not yet been sampled
                exlats = exlats[separate]
                exlons = exlons[separate]
                exvals = exvals[separate]
                exexceed = exexceed[separate]
    
                # index of max exceedance point
                imax = np.nanargmax(exexceed)
    
                # max value
                maxtp = exvals[imax]
                maxexceed = exexceed[imax]
    
                # lat/lon of max value
                maxlat = exlats[imax]
                maxlon = exlons[imax]
    
                # coordinate of max value
                maxcoord = (maxlat, maxlon)
    
                # list index of points in separate simultaneous events
                separate = []
    
                # list exceedance point lat/lons
                exceedpts = []
    
                 # check each extreme point and save index for next iteration if outside sep_dist
                for ptidx in range(len(exvals)):
                    pt = (exlats[ptidx], exlons[ptidx])
                    dist = geodesic(maxcoord, pt).km
                    if dist > sep_dist:
                        separate.append(ptidx)
                    
                # save all exceedpts within sep_dist
                # include expoints already sampled so nearby events are combined
                for ptidx in range(len(exlats_orig)):
                    pt = (exlats_orig[ptidx], exlons_orig[ptidx])
                    dist = geodesic(maxcoord, pt).km
                    if dist <= sep_dist:
                        exceedpts.append(pt)
    
                numpoints = len(exceedpts)
    
                # check previous 12 hours for time of max rain rate at point of max exceedance
                accum_1hr = []
                hour_idx = []
                for hidx in range(12):
                    accum_1hr.append(tp12[hidx][(lat==maxlat) & (lon==maxlon)][0])
                maxrate = np.nanmax(accum_1hr)
                idx_max_hour = accum_1hr.index(np.nanmax(accum_1hr))
                # fist hour within accum_period-hr window where >=1 mm of rain accumulates at point of max exceedance
                idx_start = accum_1hr.index(next(i for i in accum_1hr if i >= 1))
                # time of max precip rate (occurring within the *following* hour)
                maxtime = times12[idx_max_hour]-pd.Timedelta(1,'h')
                starttime = times12[idx_start]-pd.Timedelta(1,'h')
    
                # record data as new preliminary event
                expts_h.append(exceedpts)
                maxpt_h.append((maxlat,maxlon))
                event_time_h.append(maxtime)
                event_max_h.append(maxtp)
                event_exceedance_h.append(maxexceed)
                event_size_h.append(numpoints)
                event_state_h.append(in_us(maxlat,maxlon))
                start_time_h.append(starttime)
                accum_time_h.append(time)
                maxrate_h.append(maxrate)
                accum_array_h.append(accum)
                latver_h.append(latver)

            # exit loop if RAM exceeds 97% and will resume from this point in the month
            ram_percent = psutil.virtual_memory().percent
            if ram_percent > 97:
                times = [t for t in times if (t<=time)]
                print('\nRAM Full')
                break

    # merge hourly events with shared exceedance points within a 12-hr period
    print('\nMerging hourly events')
    total = len(maxpt_h)
    while len(maxpt_h)>0:
        # always start with first value in lists
        i=0
        
        # initialize list of event indexes that have matching points
        ematchidxs = [i]
        
        # indexes of events to check
        check_next = 500
        # can increase this value if dealing with a very large number of events in a short timespan
        irange = list(np.arange(i,i+check_next+1)[np.where(np.arange(i,i+check_next+1)<len(event_time_h))[0]])
        
        # initialize list of exceedance points for event
        expts = expts_h[i].copy()
        
        # initialize number of additional matching events
        addition = True
        
        while addition:
            addition = False
            # remove indexes already found
            irange = list(set(irange) - set(ematchidxs))
    
            if len(irange) > 0:
                # initialize i as first element in irange
                i=irange[0]
                for i2 in irange:        
                    # time difference between accumulation times
                    tdif = np.abs((accum_time_h[i] - accum_time_h[i2]) / np.timedelta64(1, 'h'))
                    tcheck = tdif < 12
            
                    # if within 12 hours check if matching exceedance points
                    if tcheck:
                        # check for common exceedance points between any match from current or previous accumulation period
                        atimeidxs = [e for e in ematchidxs if (np.abs((accum_time_h[i2] - accum_time_h[e]) / np.timedelta64(1, 'h')) <=1)]
                        expts_atime = sum([expts_h[a] for a in atimeidxs], [])
                        commonpts = common_member(expts_atime, expts_h[i2])
            
                        if len(commonpts) > 0:
                            # save index of event
                            ematchidxs.append(i2)
                            
                            # add new exceedance points
                            expts += list(set(expts_h[i2]) - set(expts))
            
                            # change i to latest match in order to continue adding events within 12 hours
                            i = i2
                            
                            addition = True
        
        # combine accumulation arrays to save only max values at each gridpoint
        matchaccums = [np.nan_to_num(accum_array_h[i]) for i in ematchidxs]
        accum_combined = np.maximum.reduce(matchaccums)
        
        # accumulation period start and end times
        period_start = np.nanmin([accum_time_h[i] for i in ematchidxs])-pd.Timedelta(12,'h')
        period_end = np.nanmax([accum_time_h[i] for i in ematchidxs])
        
        # find idx with max exceedance
        matchexceeds = [event_exceedance_h[i] for i in ematchidxs]
        emaxidx = ematchidxs[matchexceeds.index(max(matchexceeds))]
    
        event_dur = (accum_time_h[emaxidx] - start_time_h[emaxidx]) / np.timedelta64(1, 'h')
    
        maxlat,maxlon = maxpt_h[emaxidx]

        # redefine lat/lon grids for event in case of different latversion
        latver_e = int(latver_h[emaxidx])
        lat_e = np.load(latlon_dir+'stageiv_lats_'+str(latver_e)+'.npy')
        lon_e = np.load(latlon_dir+'stageiv_lons_'+str(latver_e)+'.npy')
        
        # get thresholds at point of max exceedance
        # bug fix: sometimes time of max exceedance uses a different latver if it changed during the event period
        ptidx = np.where((lat_e >= maxlat-1e-6) & (lat_e <= maxlat+1e-6) & (lon_e >= maxlon-1e-6) & (lon_e <= maxlon+1e-6))
        if len(ptidx[0])==0:
            # try each latver starting at latest
            for latver_except in range(latver,0,-1):
                lat_e = np.load(latlon_dir+'stageiv_lats_'+str(latver_except)+'.npy')
                lon_e = np.load(latlon_dir+'stageiv_lons_'+str(latver_except)+'.npy')
                ptidx = np.where((lat_e >= maxlat-1e-6) & (lat_e <= maxlat+1e-6) & (lon_e >= maxlon-1e-6) & (lon_e <= maxlon+1e-6))
                if len(ptidx[0]) > 0:
                    latver_e = latver_except
                    break
                else:
                    continue
            
        thresh = {}
        for ARI in [10,25,50,100,500,1000]:
            atlas_e = np.load(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_'+str(ARI)+'yr_regrid_'+str(int(latver_e))+'.npy')
            thresh[str(ARI)] = atlas_e[ptidx][0]
        
        # create exceedance arrays
        # mask out all points not included in this event (in case of simulataneous events)    
        exmask = (lat_e==expts[0][0]) & (lon_e==expts[0][1])
        for expt in expts[1:]:
            exmask += (lat_e==expt[0]) & (lon_e==expt[1])
        atlas_e = np.load(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_10yr_regrid_'+str(int(latver_e))+'.npy')
        exceed_array = exmask*accum_combined-atlas_e*exmask
        
        # calculate exceedance volume
        exvolume = np.nansum(exceed_array)
        
        # record data for event
        data['lat'].append(maxlat)
        data['lon'].append(maxlon)
        data['state'].append(in_us(maxlat,maxlon))
        data['event_time'].append(event_time_h[emaxidx])
        data['start_time'].append(start_time_h[emaxidx])
        data['accum_time'].append(accum_time_h[emaxidx])
        data['duration'].append(event_dur)
        data['period_start'].append(period_start)
        data['period_end'].append(period_end)
        data['exceedance'].append(event_exceedance_h[emaxidx])
        data['accumulation'].append(event_max_h[emaxidx])
        data['exceed_pts'].append(len(expts))
        data['exceed_vol'].append(exvolume)
        data['max_1hr'].append(maxrate_h[emaxidx])
        data['ari_10'].append(thresh['10'])
        data['ari_25'].append(thresh['25'])
        data['ari_50'].append(thresh['50'])
        data['ari_100'].append(thresh['100'])
        data['ari_500'].append(thresh['500'])
        data['ari_1000'].append(thresh['1000'])
        data['latversion'].append(latver_e)
        
        # append arrays to lists
        accum_arrays.append(accum_combined)
        exceed_arrays.append(exceed_array)
        exceed_points.append(expts)
        
        # remove indexes from all lists with suffix _h
        lsts = [expts_h,maxpt_h,event_time_h,event_max_h,event_exceedance_h,maxrate_h,event_size_h,event_state_h,
                start_time_h,accum_time_h,accum_array_h]
        
        for lst in lsts:
            for i in sorted(ematchidxs, reverse=True):
                del lst[i]

    # sort each list by event time
    sortedIndex = sorted(range(len(data['event_time'])), key=lambda k: data['event_time'][k])
    data_sorted = {}
    for col in cols[1:]:
        data_sorted[col] = [data[col][i] for i in sortedIndex]
    
    # for array lists
    arrs = ['accum_arrays','exceed_arrays','exceed_points']
    arrays_sorted = {}
    for arr in arrs:
        arrays_sorted[arr] = [eval(arr)[i] for i in sortedIndex]

    # append data from any events from previous period with end periods within 12 hours of the start of the period being analyzed (for merging)
    if len(df0)>0:
        df1 = df0[df0['period_end']>=times[0]-pd.Timedelta(12,'h')]
        if len(df1)>0:
            df1 = df0[df0['event']>=np.nanmin(df1['event'])]
            df0.drop(df1.index, inplace=True)
            event_id = list(df1['event'])[0]

        for col in cols[1:]:
            data_sorted[col] = list(df1[col]) + data_sorted[col]

        arrays_prev = {}
        for arr in arrs:
            arrays_prev[arr] = []
        for e in df1['event']:
            arrays_prev['accum_arrays'].append(np.load(accum_path+str(e).zfill(5)+'.npy'))
            arrays_prev['exceed_arrays'].append(np.load(exceed_path+str(e).zfill(5)+'.npy'))
            arrays_prev['exceed_points'].append([tuple(p) for p in list(np.load(exceedpt_path+str(e).zfill(5)+'.npy'))])
        for arr in arrs:
            arrays_sorted[arr] = arrays_prev[arr] + arrays_sorted[arr]
    
    # merge repeated events
    data_sorted,arrays_sorted = merge_repeats()

    # record event ids
    eids = range(event_id,event_id+len(data_sorted['lat']))
    data_sorted['event'] = eids

    print('Saving data')
    # export all arrays for the month
    for i,amap,emap,expt in zip(eids,arrays_sorted['accum_arrays'],arrays_sorted['exceed_arrays'],arrays_sorted['exceed_points']):
        np.save(accum_path+str(i).zfill(5)+'.npy',amap)
        np.save(exceed_path+str(i).zfill(5)+'.npy',emap)
        np.save(exceedpt_path+str(i).zfill(5)+'.npy',np.array(expt))
    
    # create dataframe
    df2 = pd.DataFrame(data_sorted)
    
    # move event column to front
    df2 = df2[ ['event'] + [ col for col in df2.columns if col != 'event' ] ]
    
    # merge previous
    df = pd.concat((df0,df2))

    # Convert the datetime columns to string
    df['event_time'] = df['event_time'].dt.strftime("%Y-%m-%d %H:%M:%S")
    df['accum_time'] = df['accum_time'].dt.strftime("%Y-%m-%d %H:%M:%S")
    df['period_start'] = df['period_start'].dt.strftime("%Y-%m-%d %H:%M:%S")
    df['period_end'] = df['period_end'].dt.strftime("%Y-%m-%d %H:%M:%S")
    
    # save
    df.to_csv(datadir+dbfname,index=False)
    
    print('Clearing temporary files')
    wipe_idx(path01)
    wipe_idx(path06)
    wipe_idx(path24)
    
    gc.collect()
    
    # generate and export accumulation maps
    generate_figures()

    # redefine all_times to start at the next month
    all_times = [t for t in all_times if t>times[-1]]

    current_time = timemod.strftime("%Y-%m-%d %H:%M:%S", timemod.localtime())

    elapsed_time = timemod.time() - code_start
    hours, rem = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(rem, 60)
    time_elapsed = f"Elapsed time: {int(hours):02}:{int(minutes):02}:{int(seconds):02}"
    print(time_elapsed)

    # Save to log
    with open(log_dir+'log_'+str(sep_dist)+'km_'+code_start_str+'.txt', 'a') as f:
        f.writelines('\n'+current_time)
        f.writelines('\n'+time_elapsed)
        f.writelines('\n'+'Completed up to '+times[-1].strftime("%Y-%m-%d %H:%M:%S"))
        f.writelines('\n'+'Number of events: '+str(len(df)))
        f.writelines('\n'+'Lat/lon grid version: '+str(latver)+'\n----------------------------')
        
    # exit if RAM is exceeding 80% and print/save where to resume after restarting python
    ram_percent = psutil.virtual_memory().percent
    if (ram_percent > 80) & (len(all_times)>0):
        np.save(datadir+'mem_error_log',np.array([all_times[0].strftime("%Y-%m-%d %H:%M:%S"),latver]))
        print('ERROR: Memory is almost full! Please restart Python and run code again to resume after setting "resume_from_memory_error = True"')
        mem_error = True
        break
    else:
        mem_error = False
        clear_output()

if mem_error:
    quit()
    
else:
    with open(log_dir+'log_'+str(sep_dist)+'km_'+code_start_str+'.txt', 'a') as f:
        f.writelines('\nComplete')
    print('Complete')

