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


def edyh(argv):
    # Loading mesh for run

    mesh_id    = 'meshArc4.5'
    meshpath   = '/scratch/usr/hbkoziel/mesh/'+mesh_id+'/'            # Defining path where mesh is stored
    mesh = pf.load_mesh(meshpath, usepickle=True, get3d=True)                                    # Loading mesh, stores it in mesh.****  
    #mesh = pf.fesom_mesh(meshpath, get3d=True)
    #mesh.zlevs = -mesh.zlevs                                            # Depth is made negative


    # In[7]:


    tracername = 'HEDYv2'
    first_year = int(sys.argv[1])
    last_year  = int(sys.argv[1])
    years      = np.arange(first_year,last_year+1,1)
    runid	= 'Arc12'
    months =np.linspace(0,11,12).astype(int)


    # In[8]:


    # ==============================================================================
    # Settings for netcdf file

    save_netcdf       = True      # Saves the interpolated field in netcdf file
    delete_old_netcdf = True      # If a netcdf file with the same name exists it will be deleted
    input_directory  = '/scratch/usr/hbkvsk12/hlrn3_work2/results/'+runid+'/' # Where the netcdf is saved
    output_directory  = '/scratch/usr/hbkoziel/'+runid+'/DIN_budgetv2/'


    from derivative import derivative

    bafux_2d, bafuy_2d = derivative(mesh)
    f0 = 2*7.2921e-5 * np.sin(mesh.y2/180.*np.pi)


    # In[5]:


    ncfile      = meshpath+'Arc4.5.initial.mesh.diag.nc'
    f           = Dataset(ncfile, 'r')
    NodalVol = f.variables['cluster_vol'][:]
    f.close()


    # In[12]:


    for ind in range(0,len(years)):
        depth = mesh.zlevs
        ind_max_depth = 18
        depth = depth[0:18]

        ncfile       = input_directory+runid+'.'+str(years[ind])+'.oce.mean.nc'
        f           = Dataset(ncfile, 'r')
        netcdf_name = tracername+'.'+str(years[ind])+'.monthly.nc'

        print 'year = ', years[ind]

        UVDN = np.zeros((len(months),len(mesh.x2),len(depth)-1)) * np.nan
        UDN = np.zeros((len(months),len(mesh.x2),len(depth)-1)) * np.nan
        VDN = np.zeros((len(months),len(mesh.x2),len(depth)-1)) * np.nan

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

            u_mean    = f.variables['u'][dayind,:].mean(axis=0)
            v_mean    = f.variables['v'][dayind,:].mean(axis=0)
            din_mean    = f.variables['tr01'][dayind,:].mean(axis=0)    

            uvdn = np.zeros((len(dayind1),len(mesh.x2),len(depth)-1)) * np.nan
            udn = np.zeros((len(dayind1),len(mesh.x2),len(depth)-1)) * np.nan
            vdn = np.zeros((len(dayind1),len(mesh.x2),len(depth)-1)) * np.nan
            print 'Looping over days'
            for day in range(0,len(dayind1)):
                print 'day = ', day
                u = f.variables['u'][day,:] - u_mean
                v = f.variables['v'][day,:] - v_mean
                din = f.variables['tr01'][day,:] - din_mean

                print 'Looping over depths'
                for i in range(0,len(depth)-1):
                    print 'depth = ', depth[i]
                    d_ind   = mesh.n32[:,i]-1

                    ula = np.zeros((len(mesh.x2)))* np.nan
                    vla = np.zeros((len(mesh.x2)))* np.nan
                    dinla = np.zeros((len(mesh.x2)))* np.nan 

                    # fill NaN arrays
                    ula[:] =  u[d_ind]
                    vla[:] =  v[d_ind]
                    dinla[:] =  din[d_ind]

                    # remove values with shallower depths
                    ula[ula == 0] = np.nan
                    vla[vla == 0] = np.nan
                    dinla[dinla == 0] = np.nan

                    #print 'Starting rotation'
                    (ux,vy)   = pf.vec_rotate_r2g(50, 15, -90, mesh.x2, mesh.y2, ula, vla, flag=1)

                    # derive derivative
                    uDN_elem = np.sum(ux[mesh.elem],axis=1)*np.sum(dinla[mesh.elem] * bafux_2d,axis=1)
                    vDN_elem = np.sum(vy[mesh.elem],axis=1)*np.sum(dinla[mesh.elem] * bafuy_2d,axis=1) #over triangle
                    uvDN_elem = uDN_elem + vDN_elem

                    #print 'Create node vector'
                    uvDN_node = np.zeros(len(mesh.x2))
                    uDN_node = np.zeros(len(mesh.x2))
                    vDN_node = np.zeros(len(mesh.x2))
                    ind_cnt = np.zeros(len(mesh.x2))
                    for ii in range(0,len(mesh.elem)):
                      nod_elem=mesh.elem[ii,:]
                      uDN_node[nod_elem]=uDN_node[nod_elem]+uDN_elem[ii]
                      vDN_node[nod_elem]=vDN_node[nod_elem]+vDN_elem[ii]
                      uvDN_node[nod_elem]=uvDN_node[nod_elem]+uvDN_elem[ii]
                      ind_cnt[nod_elem]=ind_cnt[nod_elem]+1.

                    # masking shallower region and 1/2 level more for noise
                    uDN_node[mesh.topo < depth[i+1]] = np.nan
                    vDN_node[mesh.topo < depth[i+1]] = np.nan
                    uvDN_node[mesh.topo < depth[i+1]] = np.nan

                    # averaging for each nodes and time convertion
                    udn[day,:,i] = 24 * 3600 * uDN_node/ind_cnt # convertion seconds to day
                    vdn[day,:,i] = 24 * 3600 * vDN_node/ind_cnt # convertion seconds to day
                    uvdn[day,:,i] = 24 * 3600 * uvDN_node/ind_cnt # convertion seconds to day

            tracershape = np.shape(UDN)

            # ==============================================================================
            # Testing if a netcdf file with the same name exists, if yes, it must be removed
            # to save a new one.

            if os.path.isfile(output_directory+netcdf_name) and delete_old_netcdf and mo==0:
              os.remove(output_directory+netcdf_name)
              print "The netcdf file "+netcdf_name+" has been deleted to make room for your file of the same name."
            elif os.path.isfile(netcdf_name) and mo==0:
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
              w_nc_fid.description = u'HEDY' 
              w_nc_fid.history     = 'Created ' + time.ctime(time.time())

              nod2d    = w_nc_fid.createDimension('nod2d', tracershape[1])  # Create dimension: number of 3d nodes
              timed	   = w_nc_fid.createDimension('time', tracershape[0]) 
              depthd    = w_nc_fid.createDimension('depth', tracershape[2]) 

              w_nc_var = w_nc_fid.createVariable('HEDY', 'f4',('time','nod2d','depth'))  # 'DIN' is name of saved variable/ 'f8' sets presicion to 64-bit floating point
              w_nc_var.setncatts({'long_name': u'HEDY',\
                                  'units': u'mmol/m2/day'})

              w_nc_var = w_nc_fid.createVariable('UEDY', 'f4',('time','nod2d','depth'))
              w_nc_var.setncatts({'long_name': u'UEDY',\
                                  'units': u'mmol/m2/day'})

              w_nc_var = w_nc_fid.createVariable('VEDY', 'f4',('time','nod2d','depth'))
              w_nc_var.setncatts({'long_name': u'VEDY',\
                                  'units': u'mmol/m2/day'})


              w_nc_fid.variables['HEDY'][mo,:,:] = np.mean(uvdn, axis = 0)   
              w_nc_fid.variables['UEDY'][mo,:,:] = np.mean(udn, axis = 0)
              w_nc_fid.variables['VEDY'][mo,:,:] = np.mean(vdn, axis = 0)  

              w_nc_fid.close()                                                     # close the new file                

              cwd = os.getcwd()
              print "New netcdf file (",netcdf_name,") has been created."
              print 'Saved month '+str(mo)+' of year '+str(years[ind])
              print "Location: "+output_directory

            elif save_netcdf and mo>0:
              w_nc_fid = Dataset(output_directory+netcdf_name, 'r+', format='NETCDF4_CLASSIC')      # Create and open new netcdf file to write to
              w_nc_fid.variables['HEDY'][mo,:,:] = np.mean(uvdn, axis = 0)   
              w_nc_fid.variables['UEDY'][mo,:,:] = np.mean(udn, axis = 0)
              w_nc_fid.variables['VEDY'][mo,:,:] = np.mean(vdn, axis = 0)  
              w_nc_fid.close()   
              print 'Saved month '+str(mo)+' of year '+str(years[ind])
            else:
              print 'You have specified not to save your field in netcdf file'
            
            t = time.localtime()
            current_time = time.strftime("%H:%M:%S", t)
            print(current_time)


edyh(sys.argv[1])

t = time.localtime()
current_time = time.strftime("%H:%M:%S", t)
print(current_time)