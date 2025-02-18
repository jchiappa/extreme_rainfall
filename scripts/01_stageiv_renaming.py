#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os


# In[ ]:


# data directories with raw (extracted) stage iv data
paths = ['E:/Research/EREs/data/stage_iv/01h/',
         'E:/Research/EREs/data/stage_iv/06h/',
         'E:/Research/EREs/data/stage_iv/24h/']


# In[ ]:


for path in paths:
    allfiles = os.listdir(path)

    # add grb2 extensions
    files = [fname for fname in allfiles if not fname.endswith('.grb2')]
    for index, file in enumerate(files):
        os.rename(path+file, path+file+'.grb2')
    
    allfiles = os.listdir(path)
    
    # start files with naming convention "st4_conus"
    files = [fname for fname in allfiles if not fname.startswith('st4_conus')]
    for index, file in enumerate(files):
        new = path+'st4_conus.'+file[-19:]
        if os.path.exists(new):
            os.remove(path+file)
        else:
            os.rename(path+file, new)

