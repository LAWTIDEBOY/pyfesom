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

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)

# In[2]:

def edyv(argv):
    # Loading mesh for run

    mesh_id    = 'meshArc4.5'
    meshpath   = '/scratch/usr/hbkoziel/mesh/'+mesh_id+'/'            # Defining path where mesh is stored
    mesh = pf.load_mesh(meshpath, usepickle=True, get3d=True)                                    # Loading mesh, stores it in mesh.****  
    #mesh = pf.fesom_mesh(meshpath, get3d=True)
    #mesh.zlevs = -mesh.zlevs                                            # Depth is made negative


    # In[3]:

    tracername = 'VEDY'
    first_year = int(sys.argv[1])
    last_year  = int(sys.argv[1])
    years      = np.arange(first_year,last_year+1,1)
    runid	= 'Arc12'
    months =np.linspace(0,11,12).astype(int)


    # In[4]:


    # ==============================================================================
    # Settings for netcdf file

    save_netcdf       = True                                            # Saves the interpolated field in netcdf file
    delete_old_netcdf = True                                            # If a netcdf file with the same name exists it will be deleted
    input_directory  = '/scratch/usr/hbkvsk12/hlrn3_work2/results/'+runid+'/' # Where the netcdf is saved
    output_directory  = '/scratch/usr/hbkoziel/'+runid+'/DIN_budgetv2/'
    plot_netcdf       = False                                         # Reads DIN from the created netcdf file, else it plots the interpolated field (should be the same)


    for ind in range(0,len(years)):
        depth = mesh.zlevs
        ind_max_depth = 18
        depth = depth[0:18]

        ncfile       = input_directory+runid+'.'+str(years[ind])+'.oce.mean.nc'
        f           = Dataset(ncfile, 'r')
        netcdf_name = tracername+'.'+str(years[ind])+'.monthly.nc'

        print 'years = ', years[ind]

        WDN = np.zeros((len(months),len(mesh.x2),len(depth)-2)) * np.nan
        print 'looping over months'
        for mo in months:
            print 'month = ',mo
            if mo ==0: 
                dayind1 = (range(0,15))
                dayind = (range(0,45))
                month='JAN'
            if mo ==1: 
                dayind1 = (range(15,30))
                dayind = (range(0,45))
                month='FEV'
            if mo ==2: 
                dayind1 = (range(30,45))
                dayind = (range(15,60))
                month='MAR'
            if mo ==3: 
                dayind1 = (range(45,60))
                dayind = (range(30,75))
                month='APR'
            elif mo ==4: 
                dayind1 = (range(60,75))
                dayind = (range(45,90))
                month='MAY'
            elif mo ==5: 
                dayind1 = (range(75,90))
                dayind = (range(60,105))
                month='JUN'
            elif mo ==6: 
                dayind1 = (range(90,105))
                dayind = (range(75,120))
                month='JUL'
            elif mo ==7:
                dayind1 = (range(105,120))
                dayind = (range(90,135))
                month='AUG'
            elif mo ==8: 
                dayind1 = (range(120,135))
                dayind = (range(105,150))
                month='SEP'
            elif mo ==9: 
                dayind1 = (range(135,150))
                dayind = (range(120,165))
                month='OCT'
            elif mo ==10: 
                dayind1 = (range(150,165))
                dayind = (range(135,167))
                month='NOV'
            elif mo ==11: 
                dayind1 = (range(165,167))
                dayind = (range(135,167))
                month='DEC'

            w_mean    = f.variables['w'][dayind,:].mean(axis=0)
            din_mean    = f.variables['tr01'][dayind,:].mean(axis=0)

            wdn = np.zeros((len(dayind1),len(mesh.x2),len(depth)-2)) * np.nan
            print 'Looping over days'
            for day in range(0,len(dayind1)):
                print 'day = ', day
                w = f.variables['w'][day,:] - w_mean
                din = f.variables['tr01'][day,:] - din_mean

                print 'Looping over nodes'
                for i in range(0,len(mesh.x2)):
                    d_ind   = mesh.n32[i,0:ind_max_depth-1]-1

                    wd = (w[d_ind[1:]]*np.diff(din[d_ind])/np.diff(depth[:ind_max_depth-1]))
                    wd[d_ind[:-1] < -99] = np.nan

                    wdn[day,i,:] = 24 * 3600 * wd

            #WDN[mo,:,:] = np.mean(wdn, axis = 0)

            tracershape = np.shape(WDN)

            # ==============================================================================
            # Testing if a netcdf file with the same name exists, if yes, it must be removed
            # to save a new one.

            if os.path.isfile(output_directory+netcdf_name) and delete_old_netcdf and mo==0:
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
            if save_netcdf and mo==0:
                import time
                w_nc_fid = Dataset(output_directory+netcdf_name, 'w', format='NETCDF4_CLASSIC')      # Create and open new netcdf file to write to
                w_nc_fid.description = u'VEDY' 
                w_nc_fid.history     = 'Created ' + time.ctime(time.time())

                nod2d    = w_nc_fid.createDimension('nod2d', tracershape[1])  # Create dimension: number of 3d nodes
                timed	   = w_nc_fid.createDimension('time', tracershape[0]) 
                depthd    = w_nc_fid.createDimension('depth', tracershape[2]) 

                w_nc_var = w_nc_fid.createVariable('VEDY', 'f4',('time','nod2d','depth'))
                w_nc_var.setncatts({'long_name': u'VEDY',\
                                      'units': u'mmol/m2/day'})

                w_nc_fid.variables['VEDY'][mo,:,:] = np.mean(wdn, axis = 0) 

                w_nc_fid.close()  # close the new file                

                cwd = os.getcwd()
                print "New netcdf file (",netcdf_name,") has been created."
                print 'Saved month '+str(mo)+' of year '+str(years[ind])
                print "Location: "+output_directory

            elif save_netcdf and mo>0:
                w_nc_fid = Dataset(output_directory+netcdf_name, 'r+', format='NETCDF4_CLASSIC')      # Create and open new netcdf file to write to
                w_nc_fid.variables['VEDY'][mo,:,:] = np.mean(wdn, axis = 0)  
                w_nc_fid.close()   
                print 'Saved month '+str(mo)+' of year '+str(years[ind])
            else:
                print 'You have specified not to save your field in netcdf file'


edyv(sys.argv[1])

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)