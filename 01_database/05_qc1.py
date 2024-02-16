#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Jason Chiappa 2024
# Adds "event duration" attribute and applies quality control conditions


# In[2]:


import os, warnings
import numpy as np
import pandas as pd

warnings.simplefilter('ignore')


# In[3]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# Main data directory
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/Publication/ERE_analysis/data/'

# Will run for entire dataset (regardless of resuming) due to short runtime

##############################################################################################


# In[4]:


df = pd.read_csv(datadir+'ere_database_prelim_qcparams.csv')

# Contert time strings to datetime
df['event_time'] = pd.to_datetime(df['event_time'], format="%Y-%m-%d %H:%M:%S")

# creates separate column for event duration in hours
df['event_duration'] = (pd.DatetimeIndex(df['accum_time']) - pd.DatetimeIndex(df['start_time'])) / np.timedelta64(1, 'h')

df


# In[5]:


# define lists from dataframe
events = df.event.to_list()
max_exceeds = df.max_exceedance.to_list()
max_accums = df.event_max.to_list()
max_sizes = df.max_size.to_list()
max_rates = df.max_rate.to_list()
surrounds = df.surround_avg.to_list()
surroundmins = df.surround_min.to_list()


# In[6]:


# function for surrounding average check
# minfraction = minimum fraction of max accumulation for surrounding gridpoints
minfraction = 0.32  # best value based on testing
def qc_suravg(event,minfraction):
    i = events.index(event)
    frac = surrounds[i]/max_accums[i]
    if frac < minfraction:
        flag = 0
    else:
        flag = 1
    return(flag)


# In[7]:


# test event
event = 11371
minfraction = 0.32
qc_suravg(event,minfraction)


# In[8]:


# function for surrounding minimum accumulation check
# minaccumfrac = minimum accumulation fraction of max exceedance accumulation for surrounding gridpoints
minaccumfrac = 0.005  # best value based on testing
def qc_surmin(event,minaccumfrac):
    i = events.index(event)
    frac = surroundmins[i]/max_accums[i]
    if frac < minaccumfrac:
        flag = 0
    else:
        flag = 1
    return(flag)


# In[9]:


# test event
event = 11371
minaccumfrac = 0.005
qc_surmin(event,minaccumfrac)


# In[10]:


# function for max accumulation rate check
# cannot exceed 104 mm/h (hail cap) - not used
# input max rainfall rate in mm/h
maxrate = 254  # 10 in/hr
def qc_rate(event,maxrate):
    i = events.index(event)
    if max_rates[i] > maxrate:
        flag = 0
    else:
        flag = 1
    return(flag)


# In[11]:


# test event
event = 11371
maxrate = 254
qc_rate(event,maxrate)


# In[12]:


# function for revmoving events with max exceedance too high for the size of the event
# input put minimum fraction of number of exceedance points to max exceedance (mm)
minfrac = 0.04  # based on testing
def qc_exceedsize(event,minfrac):
    i = events.index(event)
    if (max_sizes[i]/max_exceeds[i] < minfrac):
        flag = 0
    else:
        flag = 1
    return(flag)


# In[13]:


# test event
event = 11371
minfrac = 0.04
qc_exceedsize(event,minfrac)


# In[14]:


# lists of flags
qcsur_list = []
qcsurm_list = []
qcrate_list = []
qcexceedsize_list = []

minfraction = 0.32
minaccumfrac = 0.005
maxrate = 254
minfrac = 0.04

# test events
for event in events:
    qcsur_list.append(qc_suravg(event,minfraction))
    qcsurm_list.append(qc_surmin(event,minaccumfrac))
    qcrate_list.append(qc_rate(event,maxrate))
    qcexceedsize_list.append(qc_exceedsize(event,minfrac))


# In[15]:


df['qcsuravg'] = qcsur_list
df['qcsurmin'] = qcsurm_list
df['qcrate'] = qcrate_list
df['qcexceedsize'] = qcexceedsize_list

df


# In[16]:


# apply the four qc conditions and filter bad from dataset
df_qc = df[(df['qcsuravg']==1) & (df['qcsurmin']==1) & (df['qcrate']==1) & (df['qcexceedsize']==1)]
df_qc


# In[17]:


# list of qc flags
qcpass = df_qc['event'].to_list()
auto_qc = []
for e in events:
    if e in qcpass:
        auto_qc.append(1)
    else:
        auto_qc.append(0)


# In[18]:


# add column with qc flags
df['flag'] = auto_qc


# In[19]:


# save data with qc flags (for analysis)
df.to_csv(datadir+'ere_database_prelim_qc1_flagged.csv',index=False)


# In[20]:


# drop qc parameters and save qc-filtered dataset
df_qc1 = df_qc.drop(['qcsuravg','qcsurmin','qcrate','qcexceedsize'],axis=1)
df_qc1.to_csv(datadir+'ere_database_qc1.csv',index=False)


# In[ ]:




