#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import os, warnings
import numpy as np

warnings.simplefilter('ignore')


# In[ ]:


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
datadir = 'C:/Users/jchiappa/OneDrive - University of Oklahoma/ERE_analysis/v2/data/'

# start and end date of dataset
sdate,edate = '2002-01-01','2024-12-31'

##############################################################################################


# In[ ]:


files = os.listdir(path01)
times = pd.to_datetime(pd.Series([file[10:20] for file in files if file.endswith('.grb2')]),format='%Y%m%d%H')


# In[ ]:


alltimes = pd.Series(pd.date_range(sdate,edate, freq="h"))[:-1]


# In[ ]:


missing = []
j = 0
for i,t in enumerate(alltimes):
    if alltimes[i] == times[j]:
        j += 1
    else:
        missing.append(alltimes[i])


# In[ ]:


pd.DataFrame(missing).to_csv(datadir+'missing_01h.csv')


# In[ ]:


files = os.listdir(path06)
times = pd.to_datetime(pd.Series([file[10:20] for file in files if file.endswith('.grb2')]),format='%Y%m%d%H').to_list()


# In[ ]:


alltimes = pd.Series(pd.date_range(sdate,edate, freq="6h"))[:-1]


# In[ ]:


missing = []
for t in alltimes:
    if t not in times:
        missing.append(t)


# In[ ]:


pd.DataFrame(missing).to_csv(datadir+'missing_06h.csv')


# In[ ]:


files = os.listdir(path24)
times = pd.to_datetime(pd.Series([file[10:20] for file in files if file.endswith('.grb2')]),format='%Y%m%d%H').to_list()


# In[ ]:


alltimes = pd.Series(pd.date_range(sdate+' 12:00:00',edate+' 12:00:00', freq="24h"))[:-1]


# In[ ]:


missing = []
for t in alltimes:
    if t not in times:
        missing.append(t)


# In[ ]:


pd.DataFrame(missing).to_csv(datadir+'missing_24h.csv')

