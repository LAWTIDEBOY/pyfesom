#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from os.path import expanduser
import socket
import sys


# BELOW IS SYSTEM/PATH CONFIGURATION #
my_home = expanduser("~")
running_system = socket.gethostname()
# Configure machine's names here #
systemHLRN = 'blogin1'
systemLOCAL = 'Laurents-MacBook-Pro.local'
my_home = expanduser("~")
# Define local pass #
if running_system == systemLOCAL:
    print "computer = ", running_system
    wd_path = os.path.join(my_home,"Documents","projects", "awi-models")
elif running_system[0:1] == "b" or running_system == systemHLRN:
    print "computer = ", running_system
    wd_path = os.path.join(my_home,"awi-models")
else:
    print 'please configure your local machine : type socket.gethostname()'
print "working directory set to", wd_path
os.chdir(wd_path)
sys.path.append("codes/modules")
#########################################


# In[2]:


import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
from load_mesh_data_new import *
import numpy as np
import datetime as dt
from PiecewiseNorm import PiecewiseNorm
from netCDF4 import Dataset
import colormaps as cmaps
from matplotlib.colors import ListedColormap


# In[3]:


# PLOT CONFIGURATION #

# load colormaps #
cmap = plt.get_cmap('RdBu_r')
# figure export 
export_plot = True
# figure export definition
dpicnt=150
# choose simulation
runid='Arc12'
# choose date
first_year = 2000
last_year  = 2015
years      = np.arange(first_year,last_year+1,1)
months =np.linspace(0,11,12).astype(int)
# choose depth
get3d = True
# load the given biological tracer #
var_id1, var_id2 = 'tr06', 'tr15'

########################

# AUTOMATIC DEFINITION OF PATHS
resultpath = '/scratch/usr/hbkvsk12/hlrn3_work2/results/'+runid+'/Oldfiles/'
savepath    = wd_path+'/figures/'+runid+'/'
meshpath    = wd_path+'/data/mesh/meshArc4.5/'

# Create figure directory if it does not exist
if os.path.exists(savepath) == False:
    try:
        os.mkdir(savepath)
        print ("Creation of the directory %s successfull" % savepath)
    except OSError:
        print ("Creation of the directory %s failed" % savepath)
else:
    print ("directory %s already existing" % savepath)
    
#ncfile      = resultpath+runid+'.initial.mesh.diag.nc'
#f           = Dataset(ncfile, 'r')
#NodalArea = f.variables['cluster_area'][:]
#NodalVol    = f.variables['cluster_vol'][:]
#f.close()


# In[4]:

mesh = fesom_mesh(meshpath, get3d = get3d)
mesh.n32 = mesh.n32-1


# In[ ]:


NUMDAYS = np.empty((len(years),mesh.n2d))
DEPTHS  = np.empty((len(years),mesh.n2d))
CONC    = np.empty((len(years),mesh.n2d))

for ye in range(0,len(years)):
  numdays = np.zeros(len(mesh.x2))
  depths  = np.zeros(len(mesh.x2))
  conc    = np.zeros(len(mesh.x2))
  ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.mean.nc'
  f       = Dataset(ncfile, 'r') 

  for ind in range(106,120):# month of August
    print years[ye], ind
    chl     = f.variables[var_id1][ind,:] + f.variables[var_id2][ind,:]
    
    for i in range(0,len(mesh.x2)):
        deep_ind = mesh.n32[i,:]
        deep_ind = deep_ind[deep_ind>-10]
        loc_chl  = chl[deep_ind] # All chl in the water column
        loc_chl  = loc_chl[loc_chl>0.01]
        if (len(loc_chl)>0 and loc_chl.argmax(axis=0)>0):
          scm_i = loc_chl.argmax(axis=0)
          depths[i]  = depths[i] + mesh.zlevs[scm_i]
          numdays[i] = numdays[i] + 2.
          conc[i]    = conc[i] + loc_chl[scm_i]

  NUMDAYS[ye,:] = numdays
  DEPTHS[ye,:]  = depths/((numdays/2.))
  CONC[ye,:]    = conc/((numdays/2.))


# In[ ]:

NUMDAYS = np.nan_to_num(NUMDAYS)
CONC = np.nan_to_num(CONC)
DEPTHS = np.nan_to_num(DEPTHS)

np.savetxt(wd_path+'/Arc12/data/SCM_numdays.csv', NUMDAYS, delimiter=";")
np.savetxt(wd_path+'/Arc12/data/SCM_concentration.csv', CONC, delimiter=";")
np.savetxt(wd_path+'/Arc12/data/SCM_depth.csv', DEPTHS, delimiter=";")