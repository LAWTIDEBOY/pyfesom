#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys 
#sys.path.append("codes/modules") # add custom Vibe 's modules
#sys.path.append("../..") # add standard 's modules
# sys.path.append('/home/hbkoziel/pyfesom/codes/modules')
sys.path.append('../..') # add standard 's modules
sys.path.append('../modules')

import pyfesom as pf
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.basemap import Basemap
import numpy as np
from netCDF4 import Dataset
import os
import time


def edyv(argv):
    # Loading mesh for run

    mesh_id    = 'meshArc4.5'
    meshpath   = '/scratch/usr/hbkoziel/mesh/'+mesh_id+'/'            # Defining path where mesh is stored
    mesh = pf.load_mesh(meshpath, usepickle=True, get3d=True)                                    # Loading mesh, stores it in mesh.****  
    #mesh = pf.fesom_mesh(meshpath, get3d=True)
    #mesh.zlevs = -mesh.zlevs                                            # Depth is made negative

    tracername = 'VEDY'
    first_year = int(sys.argv[1])
    last_year  = int(sys.argv[1])
    years      = np.arange(first_year,last_year+1,1)
    runid	= 'Arc12'
    months =np.linspace(0,11,12).astype(int)

    # ==============================================================================
    # Settings for netcdf file

    save_netcdf       = True       # Saves the interpolated field in netcdf file
    delete_old_netcdf = True       # If a netcdf file with the same name exists it will be deleted
    input_directory  = '/scratch/usr/hbkvsk12/hlrn3_work2/results/'+runid+'/' # Where the netcdf is saved
    output_directory  = '/scratch/usr/hbkoziel/'+runid+'/netcdf_monthly/'
    plot_netcdf       = False       
    depth = mesh.zlevs
    depth

    # np.diff(depth)
    # mesh.n32[0,0:11]-1
    depth = mesh.zlevs
    ind_max_depth = 18
    depth[0:ind_max_depth]

    for ind in range(0,len(years)):
        ncfile       = input_directory+runid+'.'+str(years[ind])+'.oce.mean.nc'
        f           = Dataset(ncfile, 'r')
        netcdf_name = tracername+'.'+str(years[ind])+'.monthly.nc'

        print years[ind]
        # ==============================================================================
        # Loading data
        print 'Loading data'
        #w       = f.variables['w'][:,:].mean(axis=0)
        #din     = f.variables['tr01'][:,:].mean(axis=0)

        WDN = np.zeros((len(months),len(mesh.x2)))
        print 'looping over months'
        for mo in months:
            print mo
            if mo ==0: 
                dayind = (range(0,15))
                month='JAN'
            if mo ==1: 
                dayind = (range(15,30))
                month='FEV'
            if mo ==2: 
                dayind = (range(30,45))
                month='MAR'
            if mo ==3: 
                dayind = (range(45,60))
                month='APR'
            elif mo ==4: 
                dayind = (range(60,75))
                month='MAY'
            elif mo ==5: 
                dayind = (range(75,90))
                month='JUN'
            elif mo ==6: 
                dayind = (range(90,105))
                month='JUL'
            elif mo ==7:
                dayind = (range(105,120))
                month='AUG'
            elif mo ==8: 
                dayind = (range(120,135))
                month='SEP'
            elif mo ==9: 
                dayind = (range(135,150))
                month='OCT'
            elif mo ==10: 
                dayind = (range(150,165))
                month='NOV'
            elif mo ==11: 
                dayind = (range(165,167))
                month='DEC'

            w_mean = f.variables['w'][dayind,:].mean(axis=0)
            din_mean = f.variables['tr01'][dayind,:].mean(axis=0)
            count = 0
            wDN_node = np.zeros((len(mesh.x2),len(dayind)))
            for day in dayind:
                w = f.variables['w'][day,:] - w_mean
                din = f.variables['tr01'][day,:] - din_mean
                #ww = w[day,:]-w[dayind,:].mean(axis=0)
                #dd = din[day,:]-din[dayind,:].mean(axis=0)
                #wDN_node = np.zeros(len(mesh.x2))
                #print 'Looping over nodes'
                for i in range(0,len(mesh.x2)):
                    d_ind   = mesh.n32[i,0:ind_max_depth]-1
                    wDN_node[i,count] = np.mean(np.diff(w[d_ind[:]]*din[d_ind])/np.diff(depth[:ind_max_depth]))
                count = count +1

            WDN[mo,:] = np.sum(wDN_node,axis=1)* 2 * 24 * 3600 # seconds -> day, bi-daily

        tracershape = np.shape(WDN)

        # ==============================================================================
        # Testing if a netcdf file with the same name exists, if yes, it must be removed
        # to save a new one.

        if os.path.isfile(output_directory+netcdf_name) and delete_old_netcdf:
          os.remove(output_directory+netcdf_name)
          print "The netcdf file "+netcdf_name+" has been deleted to make room for your file of the same name."
        elif os.path.isfile(netcdf_name):
          statement = "The netcdf file "+netcdf_name+" already exists! It must be removed for a new one to be created. This can be done by changing your settings."
          sys.exit(statement)

        if not os.path.isdir(output_directory):
          os.makedirs(output_directory)
          print 'Directory '+output_directory+' has been created'

        # ==============================================================================
        # Creating netcdf file
        if save_netcdf:  
          import time
          w_nc_fid = Dataset(output_directory+netcdf_name, 'w', format='NETCDF4_CLASSIC')      # Create and open new netcdf file to write to
          w_nc_fid.description = u'VEDY' 
          w_nc_fid.history     = 'Created ' + time.ctime(time.time())

          nod2d    = w_nc_fid.createDimension('nod2d', mesh.n2d)               # Create dimension: number of 3d nodes
          time	   = w_nc_fid.createDimension('time', tracershape[0]) 

          w_nc_var = w_nc_fid.createVariable('VEDY', 'f4',('time','nod2d'))           # 'DIN' is name of saved variable/ 'f8' sets presicion to 64-bit floating point
          w_nc_var.setncatts({'long_name': u'VEDY',                          'units': u'mmol/m2/day'})
          w_nc_fid.variables['VEDY'][:] = WDN   


          w_nc_fid.close()                                                     # close the new file                

          cwd = os.getcwd()
          print "New netcdf file (",netcdf_name,") has been created."
          print "Location: "+output_directory
        else:
          print 'You have specified not to save your field in netcdf file'

        
edyv(sys.argv[1])