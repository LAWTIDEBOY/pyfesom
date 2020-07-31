#!/usr/bin/env python
# coding: utf-8

# In[6]:


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


# In[7]:


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


# In[8]:


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
mld_id, par_id   ='mixlay', 'PAR3D_mean'
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


# In[9]:


mesh = fesom_mesh(meshpath, get3d = get3d)
mesh.n32 = mesh.n32-1


# In[10]:


ncfile      = resultpath+runid+'.initial.mesh.diag.nc'
f           = Dataset(ncfile, 'r')
#NodalArea = f.variables['cluster_area'][:]
NodalVol3D    = f.variables['cluster_vol'][:]
f.close()


# In[ ]:


PAR  = np.zeros((len(years),len(mesh.x2)))

for ye in range(0,len(years)):
 print years[ye]
 if years[ye] < 2000:
    parml= np.zeros(len(mesh.x2))
    ind=4
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml

 
 else:
  dayind = (range(45,59))
  parml= np.zeros(len(mesh.x2))
  for ind in dayind:

    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = par[i] + np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml / len(dayind)

    
data2=PAR.mean(axis=0)
print 'Number of nans in tracer: ',np.count_nonzero(np.isnan(data2))
print 'Number of inf in tracer: ',np.count_nonzero(np.isinf(data2))
print 'Mean of surface: ',np.mean(data2[0:len(mesh.x2)])
print 'Max and min: ',np.max(data2),np.min(data2)

if export_csv == True:
    np.savetxt(outputpath+'PARmld_APR_'+str(first_year)+'_'+str(last_year)+'_trend.csv', PAR, delimiter=";")
    np.savetxt(outputpath+'PARmld_APR_'+str(first_year)+'_'+str(last_year)+'_mean.csv', data2, delimiter=";")
    print 'exporting done ...'


# In[ ]:


PAR  = np.zeros((len(years),len(mesh.x2)))

for ye in range(0,len(years)):
 print years[ye]
 if years[ye] < 2000:
    parml= np.zeros(len(mesh.x2))
    ind=5
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml

 
 else:
  dayind = (range(60,74))
  parml= np.zeros(len(mesh.x2))
  for ind in dayind:

    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = par[i] + np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml / len(dayind)

    
data2=PAR.mean(axis=0)
print 'Number of nans in tracer: ',np.count_nonzero(np.isnan(data2))
print 'Number of inf in tracer: ',np.count_nonzero(np.isinf(data2))
print 'Mean of surface: ',np.mean(data2[0:len(mesh.x2)])
print 'Max and min: ',np.max(data2),np.min(data2)

if export_csv == True:
    np.savetxt(outputpath+'PARmld_APR_'+str(first_year)+'_'+str(last_year)+'_trend.csv', PAR, delimiter=";")
    np.savetxt(outputpath+'PARmld_APR_'+str(first_year)+'_'+str(last_year)+'_mean.csv', data2, delimiter=";")
    print 'exporting done ...'


# In[ ]:


PAR  = np.zeros((len(years),len(mesh.x2)))

for ye in range(0,len(years)):
 print years[ye]
 if years[ye] < 2000:
    parml= np.zeros(len(mesh.x2))
    ind=6
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml

 
 else:
  dayind = (range(75,89))
  parml= np.zeros(len(mesh.x2))
  for ind in dayind:

    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = par[i] + np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml / len(dayind)

    
data2=PAR.mean(axis=0)
print 'Number of nans in tracer: ',np.count_nonzero(np.isnan(data2))
print 'Number of inf in tracer: ',np.count_nonzero(np.isinf(data2))
print 'Mean of surface: ',np.mean(data2[0:len(mesh.x2)])
print 'Max and min: ',np.max(data2),np.min(data2)

if export_csv == True:
    np.savetxt(outputpath+'PARmld_JUN_'+str(first_year)+'_'+str(last_year)+'_trend.csv', PAR, delimiter=";")
    np.savetxt(outputpath+'PARmld_JUN_'+str(first_year)+'_'+str(last_year)+'_mean.csv', data2, delimiter=";")
    print 'exporting done ...'


# In[ ]:


PAR  = np.zeros((len(years),len(mesh.x2)))

for ye in range(0,len(years)):
 print years[ye]
 if years[ye] < 2000:
    parml= np.zeros(len(mesh.x2))
    ind=7
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml

 
 else:
  dayind = (range(90,105))
  parml= np.zeros(len(mesh.x2))
  for ind in dayind:

    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = par[i] + np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml / len(dayind)

    
data2=PAR.mean(axis=0)
print 'Number of nans in tracer: ',np.count_nonzero(np.isnan(data2))
print 'Number of inf in tracer: ',np.count_nonzero(np.isinf(data2))
print 'Mean of surface: ',np.mean(data2[0:len(mesh.x2)])
print 'Max and min: ',np.max(data2),np.min(data2)

if export_csv == True:
    np.savetxt(outputpath+'PARmld_JUL_'+str(first_year)+'_'+str(last_year)+'_trend.csv', PAR, delimiter=";")
    np.savetxt(outputpath+'PARmld_JUL_'+str(first_year)+'_'+str(last_year)+'_mean.csv', data2, delimiter=";")
    print 'exporting done ...'


# In[ ]:


PAR  = np.zeros((len(years),len(mesh.x2)))

for ye in range(0,len(years)):
 print years[ye]
 if years[ye] < 2000:
    parml= np.zeros(len(mesh.x2))
    ind=8
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml

 
 else:
  dayind = (range(106,120))
  parml= np.zeros(len(mesh.x2))
  for ind in dayind:

    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = parml[i] + np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml / len(dayind)

    
data2=PAR.mean(axis=0)
print 'Number of nans in tracer: ',np.count_nonzero(np.isnan(data2))
print 'Number of inf in tracer: ',np.count_nonzero(np.isinf(data2))
print 'Mean of surface: ',np.mean(data2[0:len(mesh.x2)])
print 'Max and min: ',np.max(data2),np.min(data2)

if export_csv == True:
    np.savetxt(outputpath+'PARmld_AUG_'+str(first_year)+'_'+str(last_year)+'_trend.csv', PAR, delimiter=";")
    np.savetxt(outputpath+'PARmld_AUG_'+str(first_year)+'_'+str(last_year)+'_mean.csv', data2, delimiter=";")
    print 'exporting done ...'


# In[ ]:


PAR  = np.zeros((len(years),len(mesh.x2)))

for ye in range(0,len(years)):
 print years[ye]
 if years[ye] < 2000:
    parml= np.zeros(len(mesh.x2))
    ind=9
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        paml[i]    = np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml

 
 else:
  dayind = (range(121,135))
  parml= np.zeros(len(mesh.x2))
  for ind in dayind:

    ncfile  = resultpath+runid+'.'+str(years[ye])+'.bio.mean.nc'
    f       = Dataset(ncfile, 'r')    
    par     = f.variables[par_id][ind,:]
    ncfile  = resultpath+runid+'.'+str(years[ye])+'.oce.diag.nc'
    f       = Dataset(ncfile, 'r')  
    mld     = f.variables[mld_id][ind,:]
    
    for i in range(0,len(mesh.x2)):
        mld_ind = (np.abs(mld[i]-mesh.zlevs)).argmin(axis=0)
        d_ind   = mesh.n32[i,0:mld_ind]
        parml[i]    = parml[i] + np.sum(par[d_ind]* NodalVol3D[d_ind])/np.sum(NodalVol3D[d_ind])
        
    PAR[ye,:]  = parml / len(dayind)

    
data2=PAR.mean(axis=0)
print 'Number of nans in tracer: ',np.count_nonzero(np.isnan(data2))
print 'Number of inf in tracer: ',np.count_nonzero(np.isinf(data2))
print 'Mean of surface: ',np.mean(data2[0:len(mesh.x2)])
print 'Max and min: ',np.max(data2),np.min(data2)

if export_csv == True:
    np.savetxt(outputpath+'PARmld_SEP_'+str(first_year)+'_'+str(last_year)+'_trend.csv', PAR, delimiter=";")
    np.savetxt(outputpath+'PARmld_SEP_'+str(first_year)+'_'+str(last_year)+'_mean.csv', data2, delimiter=";")
    print 'exporting done ...'

