#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jason Chiappa 2025
# Organizes csv ERE database


# In[ ]:


import pandas as pd
import os, warnings

warnings.simplefilter('ignore')


# In[ ]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# Separation distance between simultaneous events
sep_dist = 250
sep_dist = 100

# Main data directory (containing ere_database csv files)
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/ERE_analysis/v2/data/'


# In[ ]:


# input database file name
dbfname = 'ere_database_classflagged_'+str(sep_dist)+'km.csv'

# output database file name (used for analysis)
outputfname = 'ere_database_'+str(sep_dist)+'km.csv'


# In[ ]:


df1 = pd.read_csv(datadir+dbfname)


# In[ ]:


# list with column names in desired order
cols = ['event','lat','lon','state','event_time','start_time','accum_time','duration','night_hours','period_start','period_end',
        'exceedance','accumulation','exceed_pts','exceed_vol','exceed_length','exceed_azimuth','max_1hr','max_pflength',
        'pflength_duration','ari_10','ari_25','ari_50','ari_100','ari_500','ari_1000','latversion','tc_flag','iso_flag','mcs_flag',
        'noct_flag','qc_flag']

df = df1[cols].convert_dtypes()


# In[ ]:


df.to_csv(datadir+outputfname,index=False)


# In[ ]:





# In[ ]:


# prep dataset for publication
df_pub = df.copy()
df_pub['event_time'] = pd.to_datetime(df_pub['event_time'], format="%Y-%m-%d %H:%M:%S")
df_pub['start_time'] = pd.to_datetime(df_pub['start_time'], format="%Y-%m-%d %H:%M:%S")
df_pub['accum_time'] = pd.to_datetime(df_pub['accum_time'], format="%Y-%m-%d %H:%M:%S")
df_pub['period_start'] = pd.to_datetime(df_pub['period_start'], format="%Y-%m-%d %H:%M:%S")
df_pub['period_end'] = pd.to_datetime(df_pub['period_end'], format="%Y-%m-%d %H:%M:%S")

syear = df_pub['event_time'].to_list()[0].year
eyear = df_pub['event_time'].to_list()[-1].year

df_pub['event_time'] = df_pub.event_time.dt.strftime('%Y%m%d%H')
df_pub['start_time'] = df_pub.start_time.dt.strftime('%Y%m%d%H')
df_pub['accum_time'] = df_pub.accum_time.dt.strftime('%Y%m%d%H')
df_pub['period_start'] = df_pub.period_start.dt.strftime('%Y%m%d%H')
df_pub['period_end'] = df_pub.period_end.dt.strftime('%Y%m%d%H')

df_pub = df_pub.drop(['duration','night_hours','latversion','exceed_azimuth'],axis=1)

outputfname = 'ere_database_'+str(syear)+'_'+str(eyear)+'_'+str(sep_dist)+'km.csv'
df_pub.to_csv(datadir+outputfname,index=False)

