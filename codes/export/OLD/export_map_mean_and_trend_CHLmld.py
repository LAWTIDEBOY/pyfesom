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
first_year = 1990
last_year  = 2015
years      = np.arange(first_year,last_year+1,1)
months =np.linspace(0,11,12).astype(int)
# choose depth
get3d = True
# load the given biological tracer #
var_id1, var_id2, mld_id = 'tr06', 'tr15', 'mixlay'
# export CSV
export_csv = True

########################

# AUTOMATIC DEFINITION OF PATHS
resultpath = '/scratch/usr/hbkvsk12/hlrn3_work2/results/'+runid+'/Oldfiles/'
savepath    = wd_path+'/figures/'+runid+'/'
meshpath    = wd_path+'/data/mesh/meshArc4.5/'
outputpath = '/scratch/usr/hbkoziel/Arc12/'

# Create figure directory if it does not exist
if os.path.exists(savepath) == False:
    try:
        os.mkdir(savepath)
        print ("Creation of the directory %s successfull" % savepath)
    except OSError:
        print ("Creation of the directory %s failed" % savepath)
else:
    print ("directory %s already existing" % savepath)


# In[4]:


mesh = fesom_mesh(meshpath, get3d = get3d)
mesh.n32 = mesh.n32-1


# In[ ]:


CHL  = np.zeros((len(years),len(mesh.x2)))

for mo in range(9,10):
    for ye in range(0,len(years)):
     print years[ye]
     if years[ye] < 2000:
        chlint= np.zeros(len(mesh.x2))
        ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.mean.nc'
        f       = Dataset(ncfile, 'r')    
        chl     = f.variables[var_id1][mo,:] + f.variables[var_id2][mo,:]
        ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
        f       = Dataset(ncfile, 'r')  
        mld     = f.variables[mld_id][mo,:]
    
        for i in range(0,len(mesh.x2)):
            # below is for weighted average
            #deep_ind   = mesh.n32[i,0:13]# integration between 0 and -190m
            #deep_ind = deep_ind[deep_ind>-10]
            #chlint[i]    = np.sum(chl[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
            # Below is for depth integration
            mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        
            if mld_ind>0:
                deep_ind = mesh.n32[i,0:mld_ind]
                loc_chl  = chl[deep_ind] # All chl in the MLD
                z_lev=mesh.z3[deep_ind]
                chlint[i] = np.trapz(loc_chl, z_lev)
            else:
                deep_ind = mesh.n32[i,mld_ind]
                loc_chl  = chl[deep_ind] # All chl in the MLD
                chlint[i] = loc_chl
      
        CHL[ye,:]  = chlint
 
     else:
      if mo ==4: 
         dayind = (range(45,59))
         month='APR'
      elif mo ==5: 
         dayind = (range(60,74))
         month='MAY'
      elif mo ==6: 
         dayind = (range(75,89))
         month='JUN'
      elif mo ==7: 
         dayind = (range(90,105))
         month='JUL'
      elif mo ==8:
         dayind = (range(106,120))
         month='AUG'
      elif mo ==9: 
         dayind = (range(121,135))
         month='SEP'
      chlint= np.zeros(len(mesh.x2))
      for ind in dayind:
          ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.mean.nc'
          f       = Dataset(ncfile, 'r')    
          chl     = f.variables[var_id1][ind,:] + f.variables[var_id2][ind,:]
          for i in range(0,len(mesh.x2)):
              mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
              if mld_ind>0:
                deep_ind = mesh.n32[i,0:mld_ind]
                loc_chl  = chl[deep_ind] # All chl in the MLD
                z_lev=mesh.z3[deep_ind]
                loc_chl=np.nan_to_num(loc_chl)
                chlint[i] = np.trapz(loc_chl, z_lev)
              else:
                deep_ind = mesh.n32[i,mld_ind]
                loc_chl  = chl[deep_ind] # All chl in the MLD
                loc_chl=np.nan_to_num(loc_chl)
                chlint[i] = loc_chl
      CHL[ye,:]  = chlint / len(dayind)

    
    data2=CHL.mean(axis=0)
    print 'Number of nans in tracer: ',np.count_nonzero(np.isnan(data2))
    print 'Number of inf in tracer: ',np.count_nonzero(np.isinf(data2))
    print 'Mean of surface: ',np.mean(data2[0:len(mesh.x2)])
    print 'Max and min: ',np.max(data2),np.min(data2)

    if export_csv == True:
        np.savetxt(outputpath+'CHLmld_'+month+'_'+str(first_year)+'_'+str(last_year)+'_trend.csv', CHL, delimiter=";")
        np.savetxt(outputpath+'CHLmld_'+month+'_'+str(first_year)+'_'+str(last_year)+'_mean.csv', data2, delimiter=";")
        print 'exporting done ...'