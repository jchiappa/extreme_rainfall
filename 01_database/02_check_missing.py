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

# st4 1h data directory
path01 = 'E:/Research/EREs/data/stage_iv/01h/'

# st4 6h data directory
path06 = 'E:/Research/EREs/data/stage_iv/06h/'

# st4 24h data directory
path24 = 'E:/Research/EREs/data/stage_iv/24h/'

# main data directory (for input and output of small files)
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/Publication/ERE_analysis/data/misc/'

# start and end date of dataset
sdate,edate = '2002-01-01','2024-01-01'

##############################################################################################


# In[3]:


files = os.listdir(path01)
times = pd.to_datetime(pd.Series([file[10:20] for file in files]),format='%Y%m%d%H')
times


# In[4]:


alltimes = pd.Series(pd.date_range(sdate,edate, freq="h"))[:-1]
alltimes


# In[5]:


missing = []
j = 0
for i,t in enumerate(alltimes):
    if alltimes[i] == times[j]:
        j += 1
    else:
        missing.append(alltimes[i])


# In[6]:


missing


# In[7]:


pd.DataFrame(missing).to_csv(datadir+'missing_01h.csv')


# In[ ]:





# In[8]:


files = os.listdir(path06)
times = pd.to_datetime(pd.Series([file[10:20] for file in files]),format='%Y%m%d%H').to_list()
times


# In[9]:


alltimes = pd.Series(pd.date_range(sdate,edate, freq="6h"))[:-1]
alltimes


# In[10]:


missing = []
for t in alltimes:
    if t not in times:
        missing.append(t)


# In[11]:


missing


# In[12]:


pd.DataFrame(missing).to_csv(datadir+'missing_06h.csv')


# In[ ]:





# In[13]:


files = os.listdir(path24)
times = pd.to_datetime(pd.Series([file[10:20] for file in files]),format='%Y%m%d%H').to_list()
times


# In[14]:


alltimes = pd.Series(pd.date_range(sdate+' 12:00:00',edate+' 12:00:00', freq="24h"))[:-1]
alltimes


# In[15]:


missing = []
for t in alltimes:
    if t not in times:
        missing.append(t)


# In[16]:


missing


# In[17]:


pd.DataFrame(missing).to_csv(datadir+'missing_24h.csv')


# In[ ]:




