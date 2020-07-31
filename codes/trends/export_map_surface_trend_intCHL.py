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


# In[ ]:


# PLOT CONFIGURATION #


# figure export definition
dpicnt=150
# choose simulation
runid='Arc12'
# choose date
first_year = 1990
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
    
ncfile      = resultpath+runid+'.initial.mesh.diag.nc'
f           = Dataset(ncfile, 'r')
NodalVol3D = f.variables['cluster_vol'][:]


# In[ ]:


mesh = fesom_mesh(meshpath, get3d = get3d)


# In[ ]:


# initlialyze
Chl2D   = np.zeros((len(years),len(mesh.x2)))
DATE = np.empty((len(years)))

# loop over years
for ind_year in range(0,len(years)):
     print years[ind_year]
     ncfile = resultpath+runid+'.'+str(years[ind_year])+'.oce.mean.nc'
     f      = Dataset(ncfile, 'r')
     chl1    = f.variables[var_id1][:,:].sum(axis=0)
     chl2    = f.variables[var_id2][:,:].sum(axis=0)
     if years[ind_year]<2000:
         Chl3D   = 30*(chl1+chl2)
     else:
         Chl3D   = 2*(chl1+chl2)
     Chl3D = Chl3D / 365
     #NPPtotal       = 365 * NPPtotal # Conversion from [mg/m2/day]   => [mg/m2/yr]
     #NPPtotal1 = np.sum(NPPtotal*NodalArea2D[:,0:mesh.n2d])/1e18 # [Pg C/year]
     for i, val in enumerate(Chl3D[0:len(mesh.x2)]):
        if i % 100000 == 0:
            print i
        ind      = mesh.n32[i,0:20]-1 # 20 = until 500m depth
        ind2     = ind[np.where(ind>-1)] # Only positive indices
        Chl2D[ind_year,i] = np.sum(Chl3D[ind2]*NodalVol3D[ind2])/np.sum(NodalVol3D[ind2])


# In[ ]:


if do_export_csv == True:
    np.savetxt("/home/hbkoziel/awi-models/data/Arc12/CHLA_annual_integrated_trend.csv", Chl2D, delimiter=";")

