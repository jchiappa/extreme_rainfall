#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Jason Chiappa 2025
# Applies quality control conditions


# In[ ]:


import os, warnings
import numpy as np
import pandas as pd

warnings.simplefilter('ignore')


# In[ ]:


##############################################################################################
#### PLEASE MODIFY AS NEEDED ####
##############################################################################################

# Separation distance between simultaneous events
sep_dist = 250
# sep_dist = 100

# Main data directory (containing ere_database_qcparams_[sep_dist]km.csv)
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/ERE_analysis/v2/data/'

# Will run for entire dataset (regardless of resuming) due to short runtime


# In[ ]:


# input database file name
dbfname = 'ere_database_qcparams_'+str(sep_dist)+'km.csv'

# output database file name
outputfname = 'ere_database_qcflagged_'+str(sep_dist)+'km.csv'

df = pd.read_csv(datadir+dbfname)


# In[ ]:


# define lists from dataframe
events = df.event.to_list()
max_exceeds = df.exceedance.to_list()
max_accums = df.accumulation.to_list()
max_sizes = df.exceed_pts.to_list()
max_rates = df.max_1hr.to_list()
surrounds = df.surround_avg.to_list()
surroundmins = df.surround_min.to_list()


# In[ ]:


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


# In[ ]:


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


# In[ ]:


# function for max accumulation rate check
# input max rainfall rate in mm/h
maxrate = 254  # 10 in/hr
def qc_rate(event,maxrate):
    i = events.index(event)
    if max_rates[i] > maxrate:
        flag = 0
    else:
        flag = 1
    return(flag)


# In[ ]:


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


# In[ ]:


# lists of flags
qcsur_list = []
qcsurm_list = []
qcrate_list = []
qcexceedsize_list = []

# minfraction = 0.32
# minaccumfrac = 0.005
# maxrate = 254
# minfrac = 0.04

# test events
for event in events:
    qcsur_list.append(qc_suravg(event,minfraction))
    qcsurm_list.append(qc_surmin(event,minaccumfrac))
    qcrate_list.append(qc_rate(event,maxrate))
    qcexceedsize_list.append(qc_exceedsize(event,minfrac))


# In[ ]:


df['qc_suravg'] = qcsur_list
df['qc_surmin'] = qcsurm_list
df['qc_rate'] = qcrate_list
df['qc_exceedsize'] = qcexceedsize_list


# In[ ]:


# apply the four qc conditions and filter bad from dataset
df_qc = df[(df['qc_suravg']==1) & (df['qc_surmin']==1) & (df['qc_rate']==1) & (df['qc_exceedsize']==1)]


# In[ ]:


# list of qc flags
qcpass = df_qc['event'].to_list()
auto_qc = []
for e in events:
    if e in qcpass:
        auto_qc.append(1)
    else:
        auto_qc.append(0)


# In[ ]:


# add column with qc flags
df['qc_flag'] = auto_qc


# In[ ]:


# save data with qc flags
df.to_csv(datadir+outputfname,index=False)

