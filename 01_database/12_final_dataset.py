#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os, warnings
import numpy as np

warnings.simplefilter('ignore')


# In[2]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# Main data directory
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/Publication/ERE_analysis/data/'

##############################################################################################


# In[3]:


df1 = pd.read_csv(datadir+'ere_database_classflagged.csv')
df1['event_time'] = pd.to_datetime(df1['event_time'], format="%Y-%m-%d %H:%M:%S")

df1


# In[4]:


df1.keys()


# In[5]:


df = pd.DataFrame({'event': df1.event,
                   'lat': df1.max_lat,
                   'lon': df1.max_lon,
                   'state': df1.state,
                   'event_time': df1.event_time,
                   'start_time': df1.start_time,
                   'accum_time': df1.accum_time,
                   'duration': df1.event_duration,
                   'night_hours': df1.night_hours,
                   'exceedance': df1.max_exceedance,
                   'accumulation': df1.event_max,
                   'exceed_pts': df1.numpts,
                   'exceed_vol': df1.ex_vol,
                   'exceed_length': df1.ex_length,
                   'exceed_azimuth': df1.ex_azimuth,
                   'max_1hr': df1.max_rate,
                   'max_pflength': df1.max_pflength,
                   'pflength_duration': df1.pflength_duration,
                   'ari_10': df1.thresh,
                   'ari_25': df1.thresh25,
                   'ari_50': df1.thresh50,
                   'ari_100': df1.thresh100,
                   'ari_500': df1.thresh500,
                   'ari_1000': df1.thresh1000,
                   'tc_flag': df1.tc_flag,
                   'iso_flag': df1.iso_flag,
                   'mcs_flag': df1.mcs_flag,
                   'noct_flag': df1.noct_flag,
                   'latversion': df1.latversion})


# In[6]:


df


# In[7]:


df.to_csv(datadir+'ere_database.csv',index=False)


# In[ ]:





# In[8]:


# prep dataset for publication
df_pub = df.copy()
df_pub['event_time'] = pd.to_datetime(df_pub['event_time'], format="%Y-%m-%d %H:%M:%S")
df_pub['start_time'] = pd.to_datetime(df_pub['start_time'], format="%Y-%m-%d %H:%M:%S")
df_pub['accum_time'] = pd.to_datetime(df_pub['accum_time'], format="%Y-%m-%d %H:%M:%S")


# In[9]:


df_pub = df_pub.drop(['duration','night_hours','state','latversion','exceed_azimuth'],axis=1)


# In[10]:


df_pub['event_time'] = df_pub.event_time.dt.strftime('%Y%m%d%H')
df_pub['start_time'] = df_pub.start_time.dt.strftime('%Y%m%d%H')
df_pub['accum_time'] = df_pub.accum_time.dt.strftime('%Y%m%d%H')


# In[11]:


df_pub = df_pub.rename(columns={'accum_time': 'end_time'})


# In[12]:


df_pub.keys()


# In[13]:


df_pub


# In[14]:


df_pub.to_csv(datadir+'ds01.csv',index=False)


# In[ ]:




