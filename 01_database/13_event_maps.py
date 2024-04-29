#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Jason Chiappa 2024
# Export image files with each event's 12-hr (right panel) and peak 1-hr (left panel) accumulation maps.
# Takes several hours and requires a lot of RAM. Memory error likely if not splitting up. Restart kernel after each portion.


# In[2]:


import os, warnings
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import shapely.geometry as sgeom
from cartopy.geodesic import Geodesic
import gc # garbage collector to clear RAM

matplotlib.rcParams['font.family'] = "Times New Roman"
matplotlib.rcParams.update({'font.size': 12})
matplotlib.rcParams['mathtext.fontset'] = "stix"

warnings.simplefilter('ignore')


# In[3]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# Enter start and end event ID number (recommend no more than 1000 for system with 32 GB of RAM)
start = 1
end = 1000

# Main data directory
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/Publication/ERE_analysis/data/'

# Accumulation map array directory
arr_dir = 'E:/Research/EREs/data/map_arrays/'

# Figure export directory
figpath = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/Publication/ERE_analysis/figs/'

# Radius from point of max exceedance to draw (km)
radius = 250

##############################################################################################

latlondir = datadir+'stageiv_latlon/'
accum_path = arr_dir+'accum_12h/'
exceed_path = arr_dir+'exceedpts/10yr/'
accum1hr_path = arr_dir+'accum_01h/'

# define export directory
figdir = figpath+'event_maps/'
os.makedirs(figdir,exist_ok=True)


# In[4]:


df = pd.read_csv(datadir+'ere_database.csv')

# Contert time strings to datetime
df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")
df['accum_time'] = pd.to_datetime(df['accum_time'], format="%Y-%m-%d %H:%M:%S")

df


# In[5]:


df.keys()


# In[6]:


events = df.event.to_list()
event_times = df.event_time.to_list()
accum_times = df.accum_time.to_list()
event_lats = df.lat.to_list()
event_lons = df.lon.to_list()
max_exceeds = df.exceedance.to_list()
max_rates = df.max_1hr.to_list()
threshs = df.ari_10.to_list()
latvers = df.latversion.to_list()
tc_flags = df.tc_flag.to_list()
iso_flags = df.iso_flag.to_list()
mcs_flags = df.mcs_flag.to_list()
noct_flags = df.noct_flag.to_list()


# In[7]:


# TEST
event = 7928
# event = 55
i = events.index(event)

lat = np.load(latlondir+'stageiv_lats_'+str(int(latvers[i]))+'.npy')
lon = np.load(latlondir+'stageiv_lons_'+str(int(latvers[i]))+'.npy')

maxlat = event_lats[i]
maxlon = event_lons[i]

thresh = threshs[i]
maxrate = max_rates[i]

tc = tc_flags[i]
iso = iso_flags[i]
mcs = mcs_flags[i]
noct = noct_flags[i]

accum1 = np.load(accum1hr_path+str(event).zfill(5)+'.npy')
# plot only 1 mm/hr pf
accum1 = np.ma.array(accum1,dtype=float,mask=accum1<1)

accum =  np.load(accum_path+str(event).zfill(5)+'.npy')
accum = np.ma.array(accum,dtype=float,mask=accum==0)

exceed =  np.load(exceed_path+str(event).zfill(5)+'.npy')

# pull out points exceeding threshold
exlats = [p[0] for p in exceed]
exlons = [p[1] for p in exceed]


# In[8]:


# create image
fig = plt.figure(figsize=(11,8),dpi=300)

# create geo-axis
ax1 = fig.add_subplot(122, projection=ccrs.PlateCarree())

# plot the data    
pcm = plt.pcolormesh(lon,lat,accum, cmap='cubehelix_r', transform=ccrs.PlateCarree(), 
                     vmin=0, vmax=thresh)

plt.scatter(exlons,exlats,marker='o',facecolors='none',edgecolors='r',s=1.5,linewidths=0.4,
            transform=ccrs.PlateCarree(),zorder=3)

plt.title(str(event).zfill(5),loc='left')
plt.title(format(accum_times[i])[:-6]+' UTC')
plt.title('{:.1f}'.format(max_exceeds[i])+' mm',loc='right')

cbar = plt.colorbar(pcm,orientation='horizontal',extend='max',pad=0.02)
cbar.set_label(label='12-hr QPE (mm)')


# create geo-axis
ax2 = fig.add_subplot(121, projection=ccrs.PlateCarree())

# plot the data    
pcm = plt.pcolormesh(lon,lat,accum1, cmap='turbo', transform=ccrs.PlateCarree(), 
                     vmin=0, vmax=maxrate)

plt.title(str(event).zfill(5),loc='left')
plt.title(format(event_times[i])[:-6]+' UTC')
plt.title('{:.1f}'.format(max_rates[i])+' mm',loc='right')

cbar = plt.colorbar(pcm,orientation='horizontal',extend='max',pad=0.02)
cbar.set_label(label='1-hr QPE (mm)')


# add geography
for ax in [ax1,ax2]:
    ax.coastlines(lw=0.5)
    ax.add_feature(cfeature.BORDERS,lw=0.5)
    ax.add_feature(cfeature.STATES,lw=0.5)

    ax.set_extent([maxlon-4.5,maxlon+4.5,maxlat-4.5,maxlat+4.5])
    
    gd = Geodesic()
    cp250 = gd.circle(lon=maxlon, lat=maxlat, radius=radius*1000)
    geom = [sgeom.Polygon(cp250)]

    ax.add_geometries(geom, crs=ccrs.PlateCarree(), facecolor='none', edgecolor='k', alpha=0.5, linewidth = 0.5)
    
    ax.scatter([maxlon],[maxlat],marker='o',facecolors='none',edgecolors='white',s=1.5,linewidths=0.6,
                transform=ccrs.PlateCarree(),zorder=3)
    
plt.subplots_adjust(wspace=0.2)

if tc==1:
    plt.text(1.1, 0.5, 'TC', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
if iso==1:
    plt.text(1.1, 0.5, 'ISO', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
if mcs==1:
    plt.text(1.1, 0.5, 'MCS', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
if noct==1:
    plt.text(1.1, 0.45, 'NOCT', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)

plt.show()


# In[9]:


if end > events[-1]:
    end = events[-1]
    
for event in range(start,end+1):
    i = events.index(event)
    
    lat = np.load(latlondir+'stageiv_lats_'+str(int(latvers[i]))+'.npy')
    lon = np.load(latlondir+'stageiv_lons_'+str(int(latvers[i]))+'.npy')

    maxlat = event_lats[i]
    maxlon = event_lons[i]

    thresh = threshs[i]
    maxrate = max_rates[i]
    
    tc = tc_flags[i]
    iso = iso_flags[i]
    mcs = mcs_flags[i]
    noct = noct_flags[i]

    accum1 = np.load(accum1hr_path+str(event).zfill(5)+'.npy')
    accum1 = np.ma.array(accum1,dtype=float,mask=accum1<1)

    accum =  np.load(accum_path+str(event).zfill(5)+'.npy')
    accum = np.ma.array(accum,dtype=float,mask=accum==0)
    
    exceed =  np.load(exceed_path+str(event).zfill(5)+'.npy')

    # pull out points exceeding threshold
    exlats = [p[0] for p in exceed]
    exlons = [p[1] for p in exceed]

    # create image
    fig = plt.figure(figsize=(11,8),dpi=300)

    # create geo-axis
    ax1 = fig.add_subplot(122, projection=ccrs.PlateCarree())

    # plot the data    
    pcm = plt.pcolormesh(lon,lat,accum, cmap='cubehelix_r', transform=ccrs.PlateCarree(), 
                         vmin=0, vmax=thresh)

    plt.scatter(exlons,exlats,marker='o',facecolors='none',edgecolors='r',s=1.5,linewidths=0.4,
                transform=ccrs.PlateCarree(),zorder=3)

    plt.title(str(event).zfill(5),loc='left')
    plt.title(format(accum_times[i])[:-6]+' UTC')
    plt.title('{:.1f}'.format(max_exceeds[i])+' mm',loc='right')

    cbar = plt.colorbar(pcm,orientation='horizontal',extend='max',pad=0.02)
    cbar.set_label(label='12-hr QPE (mm)')


    # create geo-axis
    ax2 = fig.add_subplot(121, projection=ccrs.PlateCarree())

    # plot the data    
    pcm = plt.pcolormesh(lon,lat,accum1, cmap='turbo', transform=ccrs.PlateCarree(), 
                         vmin=0, vmax=maxrate)

    plt.title(str(event).zfill(5),loc='left')
    plt.title(format(event_times[i])[:-6]+' UTC')
    plt.title('{:.1f}'.format(max_rates[i])+' mm',loc='right')

    cbar = plt.colorbar(pcm,orientation='horizontal',extend='max',pad=0.02)
    cbar.set_label(label='1-hr QPE (mm)')


    # add geography
    for ax in [ax1,ax2]:
        ax.coastlines(lw=0.5)
        ax.add_feature(cfeature.BORDERS,lw=0.5)
        ax.add_feature(cfeature.STATES,lw=0.5)

        ax.set_extent([maxlon-4.5,maxlon+4.5,maxlat-4.5,maxlat+4.5])

        gd = Geodesic()
        cp250 = gd.circle(lon=maxlon, lat=maxlat, radius=radius*1000)
        geom = [sgeom.Polygon(cp250)]

        ax.add_geometries(geom, crs=ccrs.PlateCarree(), facecolor='none', edgecolor='k', alpha=0.5, linewidth = 0.5)

        ax.scatter([maxlon],[maxlat],marker='o',facecolors='none',edgecolors='white',s=1.5,linewidths=0.6,
                    transform=ccrs.PlateCarree(),zorder=3)

    plt.subplots_adjust(wspace=0.2)
    
    if tc==1:
        plt.text(1.1, 0.5, 'TC', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
    if iso==1:
        plt.text(1.1, 0.5, 'ISO', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
    if mcs==1:
        plt.text(1.1, 0.5, 'MCS', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
    if noct==1:
        plt.text(1.1, 0.45, 'NOCT', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)

    plt.savefig(figdir+str(event).zfill(5)+'.jpg', dpi=300, bbox_inches='tight')

    fig.clf()
    plt.close(fig)
    
    # clear RAM
    gc.collect()
    
    print(str(event)+'/'+str(events[-1]),end='\r')


# In[ ]:




