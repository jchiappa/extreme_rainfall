#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Jason Chiappa 2024
# This program creates a preliminary database of extreme precipitation events from the Stage IV dataset.
# EREs defined as a single or group of 12 hr accumulation grid points exceeding the 10-yr ARI threshold from NOAA Atlas 14.
# A csv file of event attributes and numpy arrays of 12-hr precip and 1-hr precip at the peak accumulation hour are exported.
# Note: Takes around 1 hour per year of data


# In[1]:


import os, warnings
import numpy as np
import xarray as xr
from pandas import DatetimeIndex as dti
import pandas as pd
from geopy.distance import geodesic
import geopandas as gpd
from shapely.geometry import Point, Polygon
from scipy import interpolate
from math import ceil, floor
import gc # garbage collector to clear RAM

warnings.simplefilter('ignore')


# In[3]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# IMPORTANT!!!
# True if resuming (or adding additional data to append to dataset), false if starting from beginning.
# Failing to change to true will overwrite all previous data!
resume = False

# ALSO IMPORTANT!!!
# Initialize Event ID number
# 1 if running from beginning of dataset. If resuming, enter the next event ID in the preliminary dataset!
event_id = 1

# Initialize lat/lon grid version id (the grid changes slightly through the dataset)
# 1 if running for the first time. If resuming, enter previous latver entry.
latver = 1

# Enter start year and month for this run (change if resuming)
startyear = 2002
startmonth = 1
# Enter month before start (for datetime to work)
startyear1 = 2001
startmonth1 = 12

# Enter end year and month
endyear = 2023
endmonth = 12
# Last day of the month
finalday = 31

# Main data directory (for input and output of small files)
# Must include states_21basic folder
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/Publication/ERE_analysis/data/'
os.makedirs(datadir,exist_ok=True)

# st4 1h data directory
path01 = 'E:/Research/EREs/data/stage_iv/01h/'

# st4 6h data directory
path06 = 'E:/Research/EREs/data/stage_iv/06h/'

# st4 24h data directory
path24 = 'E:/Research/EREs/data/stage_iv/24h/'

# NOAA Atlas 14 data directory
atlas_dir = 'E:/Research/EREs/data/noaa_atlas_14/'

# Map array export directory (will include folders for 12-hr accumulation maps and 1-hr maps at peak accumulation hours)
arr_dir = 'E:/Research/EREs/data/map_arrays/'
os.makedirs(arr_dir,exist_ok=True)

# lat/lon bounds for domain of study (keep consistent)
west,east,south,north = -104,-65,20,50

# distance by which to separate simultaneous events (km) - i.e. radial distance from max accumulation (keep consistent)
sep_dist = 250

##############################################################################################


# In[3]:


# open shapefile for creatinge CONUS masks
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


# In[4]:


wipe_idx(path01)
wipe_idx(path06)
wipe_idx(path24)


# In[5]:


# 1h st4 files
files = [f for f in os.listdir(path01) if f.endswith('.grb2')]

# 6h st4 files
files6 = [f for f in os.listdir(path06) if f.endswith('.grb2')]

# 24h st4 files
files24 = [f for f in os.listdir(path24) if f.endswith('.grb2')]


# In[7]:


# for retrieving initial lat/lon fields
ds = xr.open_dataset(path01+'st4_conus.'+str(startyear)+str(startmonth).zfill(2)+'0101.01h.grb2')
wipe_idx(path01)

# view the dataset
ds


# In[8]:


lat = ds.latitude.data
lon = ds.longitude.data - 360


# In[9]:


# export lat/lon grids if not already done
latlondir = datadir+'stageiv_latlon/'
os.makedirs(latlondir,exist_ok=True)
if not os.path.exists(latlondir+'stageiv_lats_'+str(latver)+'.npy'):
    np.save(latlondir+'stageiv_lats_'+str(latver)+'.npy',lat)
    np.save(latlondir+'stageiv_lons_'+str(latver)+'.npy',lon)


# In[10]:


# create basic lat/lon masks for domain
latmask = (lat>=south) & (lat<=north)
lonmask = (lon>=west) & (lon<=east)
llmask = latmask & lonmask


# In[11]:


# create a location mask (within USA borders/coastlines) for the current lat/lon grid
def us_mask(lat,lon):
    global usmaskdir
    global latver
    
    os.makedirs(usmaskdir,exist_ok=True)

    binary = np.zeros(lat.shape)

    for lt in range(lat.shape[0]):
        for ln in range(lat.shape[1]):
            # check if in US
            if in_us(lat[lt,ln],lon[lt,ln]) != None:
                binary[lt,ln] = 1
                
    # location mask
    lmask = np.ma.make_mask(binary)
    
    np.save(usmaskdir+'US_mask_'+str(latver)+'.npy',lmask)
    
    return(lmask)


# In[12]:


# directory of US masks
usmaskdir = datadir+'us_masks/'

# export US mask if not already done
if not os.path.exists(usmaskdir+'US_mask_'+str(latver)+'.npy'):
    US_mask = us_mask(lat,lon)

# otherwise, load the mask file
else:
    US_mask = np.load(usmaskdir+'US_mask_'+str(latver)+'.npy')

# combine lat/lon and US masks
mask = llmask.data & US_mask
mask_i = np.invert(mask)


# In[13]:


# interpolate NOAA Atlas 14 data (12-hr precipitation frequency) onto stage iv grid
def atlas_regrid(lat,lon):
    global atlas_dir
    global latver
    
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
        atlas_nonan = np.nan_to_num(atlas_data.data,nan=1000)

        # regrid/interpolate atlas to match stage-iv
        f = interpolate.interp2d(alon,alat,atlas_nonan,kind='linear')

        atlas_regrided = np.empty_like(lat)

        for latix in range(atlas_regrided.shape[0]):
            for lonix in range(atlas_regrided.shape[1]):
                atlas_regrided[latix,lonix] = f(lon[latix,lonix],lat[latix,lonix])

        np.save(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_'+str(ARI)+'yr_regrid_'+str(latver)+'.npy',atlas_regrided)
        
        # clear RAM
        gc.collect()


# In[14]:


# load regridded 10-yr average recurrence interval (ARI) thresholds
# regrid and export if not already done
if not os.path.exists(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_10yr_regrid_'+str(latver)+'.npy'):
    atlas_regrid(lat,lon)

atlas = np.load(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_10yr_regrid_'+str(latver)+'.npy')


# In[15]:


# define array subdirectories
arr_dir_12 = arr_dir+'accum_prelim_12h/'
os.makedirs(arr_dir_12,exist_ok=True)
arr_dir_01 = arr_dir+'accum_prelim_01h/'
os.makedirs(arr_dir_01,exist_ok=True)


# In[16]:


# subroutine to retrieve 12hr precip from 1hr data files
# tp12 is the 12-hr accumulation array from the previous iteration
# tstr is the formated time string for the current iteration
def get_12hr_precip(tp12,tstr):
    global lat
    global lon
    # for first iteration, fill array with all 12 hours
    if time==times[0]:
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


# In[17]:


### use to apply 6-hr data QC ###
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


# In[18]:


### use to apply 24-hr data QC on 6-hr QC'd data###
def cross_check_24():
    global tp24sum
    # if same times as previous iteration, keep previous accumulation, else create new
    if not np.array_equal(times24,times24_old):
        # new precip array (from 24hr data)
        tp48 = np.empty((len(times24), lat.shape[0], lat.shape[1]))
        # fill array with data and sum
        for idx,time24 in enumerate(times24):
            # time in string format
            tstr = time24.strftime('%Y%m%d')+'12'
            # 24h filename
            f = 'st4_conus.'+tstr+'.24h.grb2'
            # make sure data is not missing
            if os.path.exists(path24+f):
                # open data and get 24h precip array
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


# In[19]:


# CAUTION: WILL REPLACE PREVIOUS DATA. INITIALIZATION ONLY! (make sure resume = True if resuming)
if not resume:
    # list of times with max hourly total (in max gridpoint) per event
    event_time = []

    # list of lats/lons/state of max accum_period-hr total per event
    event_lat = []
    event_lon = []
    event_state = []

    # list of values of max accum_period-hr totals per event
    event_max = []

    # list of maximum 1-hr rainfall at point of max exceedance
    event_max_rate = []

    # exceedance threshold at point
    event_thresh = []

    # amount exceeded
    event_exceedance = []

    # end time of accum_period-hr accumulation interval with max value
    accum_time = []
    # start time of precip within accum_period-hour window at point of max exceedance
    start_time = []
    # list max accum_period-hr total # of gridpoints exceeding threshold per event
    event_size = []

    # list of maps for each event
    figs_accum = []
    figs_accum1hr = []

    # list of lat/lon grid version IDs
    latlonvers = []

    eids0 = range(0,len(event_time))

    df = pd.DataFrame({'event': eids0, 'event_time': event_time, 'start_time': start_time, 'accum_time': accum_time, 
                           'event_max': event_max, 'thresh': event_thresh, 'max_exceedance': event_exceedance, 'max_size': event_size, 
                           'max_lat': event_lat, 'max_lon': event_lon, 'max_rate': event_max_rate, 'state': event_state,
                           'latversion': latlonvers})

    df.to_csv(datadir+'ere_database_prelim.csv',index=False)


# In[21]:


# create series of datetimes for each month to collect data from
starttime = str(startyear)+'-'+str(startmonth).zfill(2)+'-01 00:00:00'
starttime1 = str(startyear1)+'-'+str(startmonth1).zfill(2)+'-01 00:00:00'
endtime = str(endyear)+'-'+str(endmonth).zfill(2)+'-'+str(finalday)+' 23:00:00'

monthstarttimes = (pd.Series(pd.date_range(starttime1,endtime, freq='m'))+pd.Timedelta(24,unit='h'))[:-1]
monthendtimes = pd.Series(pd.date_range(starttime,endtime, freq='m'))+pd.Timedelta(23,unit='h')

# initialize cross-check file times
times06_old = 0
times24_old = 0
    
# collect events one month at a time
for starttime,endtime in zip(monthstarttimes,monthendtimes):
    # print out year/month currently running
    print(str(starttime.month)+'/'+str(starttime.year))

    # dataframe up to previous month
    df0 = pd.read_csv(datadir+'ere_database_prelim.csv')
    if resume:
        df0['event_time'] = pd.to_datetime(df0['event_time'], format="%Y-%m-%d %H:%M:%S")
        # save parameters for events before run period
        df0 = df0[df0.event_time < starttime]

    # times to sum previous 12 hours of data
    times = pd.Series(pd.date_range(starttime,endtime, freq='h'))

    ### Define empty lists to store data
    # list of times with max hourly total (in max gridpoint) per event
    event_time = []

    # list of lats/lons/state of max accum_period-hr total per event
    event_lat = []
    event_lon = []
    event_state = []

    # list of values of max accum_period-hr totals per event
    event_max = []

    # list of maximum 1-hr rainfall at point of max exceedance
    event_max_rate = []

    # exceedance threshold at point
    event_thresh = []

    # amount exceeded
    event_exceedance = []

    # end time of accum_period-hr accumulation interval with max value
    accum_time = []
    # start time of precip within accum_period-hour window at point of max exceedance
    start_time = []
    # list max accum_period-hr total # of gridpoints exceeding threshold per event
    event_size = []

    # list of maps for each event
    figs_accum = []
    figs_accum1hr = []

    # list of lat/lon grid version IDs
    latlonvers = []

    # 12hr precip array
    tp12 = np.empty((12, lat.shape[0], lat.shape[1]))

    # loop for each hour in the month
    for index,time in enumerate(times):
        # show progress percentage for the month
        print(str(time)+' '+"{:.2f}".format((index+1)/len(times)*100)+'%',end='\r')
        
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
            if not os.path.exists(latlondir+'stageiv_lats_'+str(latver)+'.npy'):
                np.save(latlondir+'stageiv_lats_'+str(latver)+'.npy',lat)
                np.save(latlondir+'stageiv_lons_'+str(latver)+'.npy',lon)

            latmask = (lat>=south) & (lat<=north)
            lonmask = (lon>=west) & (lon<=east)
            llmask = latmask & lonmask

            # redo usmask
            US_mask = us_mask(lat,lon)
            
            mask = llmask.data & US_mask
            mask_i = np.invert(mask)
            
            # regrid atlas data
            atlas_regrid(lat,lon)
            atlas = np.load(atlas_dir+'NOAA_Atlas_14_CONUS_pf_012_hr_10yr_regrid_'+str(latver)+'.npy')
            
        # perform cross-checks on data to take advantage of manual QC on 6 and 24 hr data
        accum_qc06 = cross_check_06()
        accum = cross_check_24()        
        
        # save times for next iteration to avoid unnecessary repeated opening of files
        times06_old = times06.copy()
        times24_old = times24.copy()

        # apply mask
        accum[mask_i] = np.nan

        # subtract threshold map
        exceed = accum-atlas

        # if any points exceed threshold
        if np.nanmax(exceed) > 0:
            # pull out points exceeding threshold
            extreme_mask = exceed>0
            exlats = lat[extreme_mask]
            exlons = lon[extreme_mask]
            exvals = accum[extreme_mask]
            exthresh = atlas[extreme_mask]
            exexceed = exceed[extreme_mask]

            # check for simultaneous events within sep_dist km
            # list index of points in separate simultaneous events (all points for first iteration)
            separate = list(range(len(exvals)))

            # while points not sampled yet, count additional events
            while len(separate) > 0:
                # index points that have not yet been sampled
                exlats = exlats[separate]
                exlons = exlons[separate]
                exvals = exvals[separate]
                exthresh = exthresh[separate]
                exexceed = exexceed[separate]

                # index of max exceedance point
                imax = np.nanargmax(exexceed)

                # max value
                maxtp = exvals[imax]
                maxthresh = exthresh[imax]
                maxexceed = exexceed[imax]

                # lat/lon of max value
                maxlat = exlats[imax]
                maxlon = exlons[imax]

                # coordinate of max value
                maxcoord = (maxlat, maxlon)

                # list index of points in separate simultaneous events
                separate = []

                # list index of points in event of interest for this sample
                eventpts = []

                # check each extreme point and save index if within sep_dist
                for ptidx in range(len(exvals)):
                    pt = (exlats[ptidx], exlons[ptidx])
                    dist = geodesic(maxcoord, pt).km
                    if dist > sep_dist:
                        separate.append(ptidx)
                    else:
                        eventpts.append(ptidx)

                numpoints = len(eventpts)

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
                # array of accumulation at max hour
                maxaccumarr = tp12[idx_max_hour]

                # index of current record
                index_current = len(event_time)

                replace_index = []

                # check for other events within accum_period hours and sep_dist and save only one with max exceedance/size
                if index_current > 0:
                    # check previous 10 data points
                    for index in list(range(index_current-10,index_current)):
                        if index < 0:
                            index = 0
                        # time difference from previous record (hours)
                        tdelta = (maxtime-event_time[index]) / np.timedelta64(1, 'h')
                        # no 2 events in same area within 12 hours
                        if tdelta <= 12:
                            # distance (km)
                            pt1 = (event_lat[index], event_lon[index])
                            dist = geodesic(maxcoord, pt1).km
                            if dist <= sep_dist:
                                # difference in number of exceedance points
                                dsize = numpoints-event_size[index]
                                # difference in max exceedance value
                                dmax = maxexceed-event_exceedance[index]

                                # to replace previous data point
                                if (dmax>0) | ((dmax >= 0) & (dsize > 0)):
                                    # so that only one can be replaced
                                    if sum(replace_index) == 0:
                                        # in case it is the 0 index to be replaced
                                        if index == 0:
                                            replace_index.append(-1)
                                        else:
                                            replace_index.append(index)

                                # to replace size only
                                elif dsize > 0:
                                    # so that only one can be replaced
                                    if sum(replace_index) == 0:
                                        # in case it is the 0 index to be replaced
                                        if index == 0:
                                            replace_index.append(-1e6)
                                        else:
                                            replace_index.append(index*1e6)

                                # to record no data    
                                else:
                                    replace_index.append(0)

                # conditions to record new data
                if len(replace_index) == 0:
                    event_time.append(maxtime)
                    event_lat.append(maxlat)
                    event_lon.append(maxlon)
                    event_max.append(maxtp)
                    event_thresh.append(maxthresh)
                    event_exceedance.append(maxexceed)
                    event_size.append(numpoints)
                    event_state.append(in_us(maxlat,maxlon))
                    start_time.append(starttime)
                    accum_time.append(time)
                    event_max_rate.append(maxrate)
                    latlonvers.append(latver)

                    # append map arrays to list for plotting
                    figs_accum.append(accum)
                    figs_accum1hr.append(maxaccumarr)

                # replace size only
                elif np.abs(sum(replace_index)) >= 1e6:
                    # index to replace
                    replace = int(sum(replace_index)/1e6)
                    # if it is 0 index
                    if replace < 0:
                        replace = 0
                    event_size[replace] = numpoints

                # conditions to replace previous data point
                elif sum(replace_index) != 0:
                    # index to replace
                    replace = sum(replace_index)
                    if replace < 0:
                        replace = 0
                    event_time[replace] = maxtime
                    event_lat[replace] = maxlat
                    event_lon[replace] = maxlon
                    event_max[replace] = maxtp
                    event_thresh[replace] = maxthresh
                    event_exceedance[replace] = maxexceed
                    event_size[replace] = numpoints
                    event_state[replace] = in_us(maxlat,maxlon)
                    start_time[replace] = starttime
                    accum_time[replace] = time
                    event_max_rate[replace] = maxrate
                    latlonvers[replace] = latver

                    figs_accum[replace] = accum
                    figs_accum1hr[replace] = maxaccumarr

                # conditions to not record data
                # if len(replace_index) > 0 and sum == 0

        # clear RAM
        gc.collect()

    eids = range(event_id,event_id+len(event_time))

    for i,amap,a1map in zip(eids,figs_accum,figs_accum1hr):
        np.save(arr_dir_12+str(i).zfill(5)+'.npy',amap)
        np.save(arr_dir_01+str(i).zfill(5)+'.npy',a1map)

    # clear RAM
    del(figs_accum)
    gc.collect()
    del(figs_accum1hr)
    gc.collect()

    df1 = pd.DataFrame({'event': eids, 'event_time': event_time, 'start_time': start_time, 'accum_time': accum_time, 
                       'event_max': event_max, 'thresh': event_thresh, 'max_exceedance': event_exceedance, 'max_size': event_size, 
                       'max_lat': event_lat, 'max_lon': event_lon, 'max_rate': event_max_rate, 'state': event_state,
                       'latversion': latlonvers})

    df = pd.concat((df0,df1))

    df.to_csv(datadir+'ere_database_prelim.csv',index=False)

    event_id = event_id+len(event_time)
    # in case of error, print event id and latver to resume from
    print('Next Event ID: '+str(event_id),' Lat/lon grid version: '+str(latver))

    wipe_idx(path01)
    wipe_idx(path06)
    wipe_idx(path24)


# In[22]:


wipe_idx(path01)
wipe_idx(path06)
wipe_idx(path24)


# In[ ]:




