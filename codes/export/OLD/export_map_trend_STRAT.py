#!/usr/bin/env python
# coding: utf-8

# In[39]:


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
elif running_system[0] == "b" or running_system == systemHLRN:
    print "computer = ", running_system
    wd_path = os.path.join(my_home,"awi-models")
else:
    print 'please configure your local machine : type socket.gethostname()'
print "working directory set to", wd_path
os.chdir(wd_path)
sys.path.append("codes/modules")
#########################################
# APR = 


# In[40]:


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


mesh = fesom_mesh(meshpath, get3d = get3d)


def dens0(T,S):
    T68 = T * 1.00024
    a0,a1,a2,a3,a4,a5 = 999.842594, 6.793952e-2, -9.095290e-3, 1.001685e-4, -1.120083e-6, 6.536332e-9
    sw_smow = a0 + (a1 + (a2 + (a3 + (a4 + a5*T68)*T68)*T68)*T68)*T68
    #     UNESCO 1983 eqn(13) p17.
    
    b0, b1, b2, b3, b4 =  8.24493e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9
    c0, c1, c2 = -5.72466e-3, 1.0227e-4, -1.6546e-6
    d0 = 4.8314e-4
    dens = sw_smow + (b0 + (b1 + (b2 + (b3 + b4*T68)*T68)*T68)*T68)*S + (c0 + (c1 + c2*T68)*T68)*S*np.sqrt(S) + d0*S**2
    return dens


# initlialyze
STRAT = np.empty((len(years),mesh.n2d))

# loop over years
for ind in range(0,len(years)):
     print years[ind]
     ncfile = resultpath+runid+'.'+str(years[ind])+'.oce.mean.nc'
     f      = Dataset(ncfile, 'r')
     if years[ind]<2000:
         salt    = f.variables['salt'][4,:]
         temp    = f.variables['temp'][4,:]
     else:
         salt    = f.variables['salt'][45:59,:].mean(axis=0)#75:135
         temp    = f.variables['temp'][45:59,:].mean(axis=0)
            
     dens = dens0(temp,salt)
     del salt
     del temp
     dens_0m        = fesom2depth(0, dens, mesh)
     dens_100m      = fesom2depth(100, dens, mesh)
    
     STRAT[ind,:]=dens_100m-dens_0m
     del dens_0m
     del dens_100m
        
if export_csv == True:
    np.savetxt(wd_path+'/data/Arc12/STRAT_APR_1990_2015_trend.csv', STRAT, delimiter=";")
print 'exporting April done...'
    
# initlialyze
STRAT = np.empty((len(years),mesh.n2d))

# loop over years
for ind in range(0,len(years)):
     print years[ind]
     ncfile = resultpath+runid+'.'+str(years[ind])+'.oce.mean.nc'
     f      = Dataset(ncfile, 'r')
     if years[ind]<2000:
         salt    = f.variables['salt'][5,:]
         temp    = f.variables['temp'][5,:]
     else:
         salt    = f.variables['salt'][60:74,:].mean(axis=0)#75:135
         temp    = f.variables['temp'][60:74,:].mean(axis=0)
            
     dens = dens0(temp,salt)
     del salt
     del temp
     dens_0m        = fesom2depth(0, dens, mesh)
     dens_100m      = fesom2depth(100, dens, mesh)
    
     STRAT[ind,:]=dens_100m-dens_0m
     del dens_0m
     del dens_100m
        
if export_csv == True:
    np.savetxt(wd_path+'/data/Arc12/STRAT_MAY_1990_2015_trend.csv', STRAT, delimiter=";")
print 'exporting May done...'
    
# initlialyze
STRAT = np.empty((len(years),mesh.n2d))

# loop over years
for ind in range(0,len(years)):
     print years[ind]
     ncfile = resultpath+runid+'.'+str(years[ind])+'.oce.mean.nc'
     f      = Dataset(ncfile, 'r')
     if years[ind]<2000:
         salt    = f.variables['salt'][6,:]
         temp    = f.variables['temp'][6,:]
     else:
         salt    = f.variables['salt'][75:89,:].mean(axis=0)#75:135
         temp    = f.variables['temp'][75:89,:].mean(axis=0)
            
     dens = dens0(temp,salt)
     del salt
     del temp
     dens_0m        = fesom2depth(0, dens, mesh)
     dens_100m      = fesom2depth(100, dens, mesh)
    
     STRAT[ind,:]=dens_100m-dens_0m
     del dens_0m
     del dens_100m
        
if export_csv == True:
    np.savetxt(wd_path+'/data/Arc12/STRAT_JUN_1990_2015_trend.csv', STRAT, delimiter=";")
print 'exporting June done...'

# initlialyze
STRAT = np.empty((len(years),mesh.n2d))

# loop over years
for ind in range(0,len(years)):
     print years[ind]
     ncfile = resultpath+runid+'.'+str(years[ind])+'.oce.mean.nc'
     f      = Dataset(ncfile, 'r')
     if years[ind]<2000:
         salt    = f.variables['salt'][7,:]
         temp    = f.variables['temp'][7,:]
     else:
         salt    = f.variables['salt'][90:105,:].mean(axis=0)#75:135
         temp    = f.variables['temp'][90:105,:].mean(axis=0)
            
     dens = dens0(temp,salt)
     del salt
     del temp
     dens_0m        = fesom2depth(0, dens, mesh)
     dens_100m      = fesom2depth(100, dens, mesh)
    
     STRAT[ind,:]=dens_100m-dens_0m
     del dens_0m
     del dens_100m
        
if export_csv == True:
    np.savetxt(wd_path+'/data/Arc12/STRAT_JUL_1990_2015_trend.csv', STRAT, delimiter=";")
print 'exporting July done...'
    
# initlialyze
STRAT = np.empty((len(years),mesh.n2d))

# loop over years
for ind in range(0,len(years)):
     print years[ind]
     ncfile = resultpath+runid+'.'+str(years[ind])+'.oce.mean.nc'
     f      = Dataset(ncfile, 'r')
     if years[ind]<2000:
         salt    = f.variables['salt'][8,:]
         temp    = f.variables['temp'][8,:]
     else:
         salt    = f.variables['salt'][106:120,:].mean(axis=0)#75:135
         temp    = f.variables['temp'][106:120,:].mean(axis=0)
            
     dens = dens0(temp,salt)
     del salt
     del temp
     dens_0m        = fesom2depth(0, dens, mesh)
     dens_100m      = fesom2depth(100, dens, mesh)
    
     STRAT[ind,:]=dens_100m-dens_0m
     del dens_0m
     del dens_100m
        
if export_csv == True:
    np.savetxt(wd_path+'/data/Arc12/STRAT_AUG_1990_2015_trend.csv', STRAT, delimiter=";")
print 'exporting August done...'
    
# initlialyze
STRAT = np.empty((len(years),mesh.n2d))

# loop over years
for ind in range(0,len(years)):
     print years[ind]
     ncfile = resultpath+runid+'.'+str(years[ind])+'.oce.mean.nc'
     f      = Dataset(ncfile, 'r')
     if years[ind]<2000:
         salt    = f.variables['salt'][9,:]
         temp    = f.variables['temp'][9,:]
     else:
         salt    = f.variables['salt'][121:135,:].mean(axis=0)#75:135
         temp    = f.variables['temp'][121:135,:].mean(axis=0)
            
     dens = dens0(temp,salt)
     del salt
     del temp
     dens_0m        = fesom2depth(0, dens, mesh)
     dens_100m      = fesom2depth(100, dens, mesh)
    
     STRAT[ind,:]=dens_100m-dens_0m
     del dens_0m
     del dens_100m
        
if export_csv == True:
    np.savetxt(wd_path+'/data/Arc12/STRAT_SEP_1990_2015_trend.csv', STRAT, delimiter=";")
print 'exporting September done...'