#!/usr/bin/env python
# coding: utf-8

# In[34]:


import os
from os.path import expanduser
import socket
import sys
import glob

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
    wd_path = os.path.join(my_home,"Documents","projects", "pyfesom")
elif running_system[0:1] == "b" or running_system == systemHLRN:
    print "computer = ", running_system
    wd_path = os.path.join(my_home,"pyfesom")
else:
    print 'please configure your local machine : type socket.gethostname()'
print "working directory set to", wd_path
os.chdir(wd_path)
sys.path.append("codes/modules") # add custom Vibe 's modules
#sys.path.append("pyfesom") # add standard 's modules
#########################################


# In[35]:


import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from pylab import *
import pyfesom as pf
#from load_mesh_data_new import *
import numpy as np
import datetime as dt
from PiecewiseNorm import PiecewiseNorm
from netCDF4 import Dataset
import cmocean as cmaps
import pandas as pd
#from matplotlib.colors import ListedColormap


# In[36]:


# PLOT CONFIGURATION #

# load colormaps #
cmap = cmaps.cm.balance
# figure export 
export_plot = True
# figure export definition
dpicnt=150
# choose simulation
runid='Arc12'
# choose date
first_year = 1980
last_year  = 2015
years      = np.arange(first_year,last_year+1,1)
months =np.linspace(0,11,12).astype(int)
# choose depth
get3d = False
# load the given biological tracer #
var_id= 'area'

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


# In[37]:


mesh = pf.load_mesh(meshpath, usepickle=True,get3d=get3d)
#mesh = fesom_mesh(meshpath, get3d = get3d)
#mesh = load_mesh(meshpath)


# In[38]:


mesh


# In[39]:


ncfile      = meshpath+'Arc4.5.initial.mesh.diag.nc'
f           = Dataset(ncfile, 'r')
NodalAreaModel = f.variables['cluster_area'][:]


# In[40]:


#X = genfromtxt(outputpath+'ICE_APR_1990_2015_trend.csv', delimiter=';')
#fig, axs = plt.subplots(3,2, sharex=True, sharey=True)

for mo in range(8,9):
    if mo ==3: 
        dayind = (range(45,59))
        month='APR'
    elif mo ==4: 
        dayind = (range(60,74))
        month='MAY'
    elif mo ==5: 
        dayind = (range(75,89))
        month='JUN'
    elif mo ==6: 
        dayind = (range(90,105))
        month='JUL'
    elif mo ==7:
        dayind = (range(106,120))
        month='AUG'
    elif mo ==8: 
        dayind = (range(121,135))
        month='SEP'
    
    # initlialyze
    ICECON = np.zeros((len(years),len(mesh.x2)))
    ICEAREA = np.zeros((len(years),len(mesh.x2)))

    # loop over years
    for ind in range(0,len(years)):
         print years[ind]
         ncfile = resultpath+runid+'.'+str(years[ind])+'.ice.mean.nc'
         #print ncfile
         f      = Dataset(ncfile, 'r')
         if years[ind]<2000:
             ice    = f.variables[var_id][mo,0:mesh.n2d]
         else:
             ice    = f.variables[var_id][dayind,0:mesh.n2d].mean(axis=0)
             #ice    = ice / (len(dayind))
         arcsurf_ind = np.nonzero((mesh.y2>65.) & (ice>0.1))
         ICEAREA[ind,:]=np.sum(NodalAreaModel[arcsurf_ind])
         ICECON[ind,:] =np.mean(ice[mesh.y2>66.])

    data2 = ICECON.mean(axis=0)
    print 'Number of nans in tracer: ',np.count_nonzero(np.isnan(data2))
    print 'Number of inf in tracer: ',np.count_nonzero(np.isinf(data2))
    print 'Mean of surface: ',np.mean(data2[0:len(mesh.x2)])
    print 'Max and min: ',np.max(data2),np.min(data2)
    
    if mo < 10:
        BSOfile = wd_path+'/data/NSIDC/N_0'+str(mo)+'_extent_v3.0.csv'
        print BSOfile
        df = pd.read_csv(BSOfile, header=0, sep = ',',na_values=['-9999'], decimal=".",skipinitialspace=True)
    else:
        BSOfile = wd_path+'/data/NSIDC/N_'+str(mo)+'_extent_v3.0.csv'
        print BSOfile
        df = pd.read_csv(BSOfile, header=0, sep = ',',na_values=['-9999'], decimal=".",skipinitialspace=True)
    icearea=ICEAREA.mean(axis=1)/1e6
    plt.plot(years,icearea,'.-')
    plt.plot(df.year,df.extent*1e6,'.-r')
    plt.grid(True)
    plt.xlabel('time (years)')
    plt.ylabel('Sea-ice extent ($km^2$)')
    plt.title(month)
    legend(['FESOM','NSIDC'])
#df.head() 


# In[ ]:





# In[41]:



#plt.savefig(savepath+'sertemp/SEAICEAREA_trend_'+str(first_year)+'_'+str(last_year)+'.png', dpi = dpicnt, bbox_inches='tight') 
     

