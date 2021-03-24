#!/usr/bin/env python
# coding: utf-8

import sys 
sys.path.append('../..') # add standard 's modules
sys.path.append('../modules')
import pyfesom as pf
import numpy as np
from netCDF4 import Dataset
import os
import time


def edyhv2(argv):
    # Loading mesh for run

    mesh_id    = 'meshArc4.5'
    meshpath   = '/scratch/usr/hbkoziel/mesh/'+mesh_id+'/'            # Defining path where mesh is stored
    mesh = pf.load_mesh(meshpath, usepickle=True, get3d=True)                                    # Loading mesh, stores it in mesh.****  
    #mesh = pf.fesom_mesh(meshpath, get3d=True)
    #mesh.zlevs = -mesh.zlevs                                            # Depth is made negative

    tracername = 'HEDYv2'
    first_year = int(sys.argv[1])
    last_year  = int(sys.argv[1])
    print first_year
    print last_year
    print 'calculating time vector...'
    #x = range(1985,1985+1,1)
    #print x
    years = np.arange(first_year,last_year+1,1)
    print years
    runid	= 'Arc12'
    print runid
    months = np.linspace(0,11,12).astype(int)
    print months

    print 'print flags'
    save_netcdf       = True           # Saves the interpolated field in netcdf file
    delete_old_netcdf = True           # If a netcdf file with the same name exists it will be deleted
    
    print 'print directories'
    input_directory  = '/scratch/usr/hbkvsk12/hlrn3_work2/results/'+runid+'/' # Where the netcdf is saved
    output_directory  = '/scratch/usr/hbkoziel/'+runid+'/netcdf_monthly/'
    
    print 'calculating derivative'
    from derivative import derivative

    bafux_2d, bafuy_2d = derivative(mesh)
    f0 = 2*7.2921e-5 * np.sin(mesh.y2/180.*np.pi)


    ncfile      = meshpath+'Arc4.5.initial.mesh.diag.nc'
    f           = Dataset(ncfile, 'r')
    NodalVol = f.variables['cluster_vol'][:]
    f.close()
    print 'diagnostic loaded'
    
    depth = mesh.zlevs
    ind_max_depth = 18
    depth[0:ind_max_depth]


    for ind in range(0,len(years)):
        ncfile       = input_directory+runid+'.'+str(years[ind])+'.oce.mean.nc'
        print(ncfile)
        f           = Dataset(ncfile, 'r')
        netcdf_name = tracername+'.'+str(years[ind])+'.monthly.nc'
        print(netcdf_name)
        print years[ind]
        print 'Loading data'
        UVDN = np.zeros((len(months),len(mesh.x2)))
        UDN = np.zeros((len(months),len(mesh.x2)))
        VDN = np.zeros((len(months),len(mesh.x2)))
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
            
            print dayind
            ula = np.zeros((len(mesh.x2),len(dayind)))
            vla = np.zeros((len(mesh.x2),len(dayind)))
            dinla = np.zeros((len(mesh.x2),len(dayind)))
            print 'test'
            u_mean = f.variables['u'][:,:]
            print np.shape(u_mean)
            print 'derive monthly means'
            u_mean = f.variables['u'][dayind,:].mean(axis=0)
            v_mean = f.variables['v'][dayind,:].mean(axis=0)
            din_mean = f.variables['tr01'][dayind,:].mean(axis=0)
            print 'loop over days'
            count = 0
            for day in dayind:
                print day
                u = f.variables['u'][day,:] - u_mean
                v = f.variables['v'][day,:] - v_mean
                din = f.variables['tr01'][day,:] - din_mean
                print 'loop over nodes'
                for i in range(0,len(mesh.x2)):
                    d_ind   = mesh.n32[i,0:ind_max_depth]-1
                    ula[i,count] = np.sum(u[d_ind]* NodalVol[d_ind])/np.sum(NodalVol[d_ind])
                    vla[i,count] = np.sum(v[d_ind]* NodalVol[d_ind])/np.sum(NodalVol[d_ind])
                    dinla[i,count] = np.sum(din[d_ind]* NodalVol[d_ind])/np.sum(NodalVol[d_ind])
                count = count +1

            ula = np.sum(ula,axis=1)*2* 24 * 3600 # seconds -> day, bi-daily
            vla = np.sum(vla,axis=1)*2* 24 * 3600 # seconds -> day, bi-daily
            print 'Starting rotation'
            (u,v)   = pf.vec_rotate_r2g(50, 15, -90, mesh.x2, mesh.y2, ula, vla, flag=1)
            uDN_elem = np.sum(u[mesh.elem],axis=1)*np.sum(din[mesh.elem] * bafux_2d,axis=1)
            vDN_elem = np.sum(v[mesh.elem],axis=1)*np.sum(din[mesh.elem] * bafuy_2d,axis=1) #over triangle
            uvDN_elem = uDN_elem + vDN_elem

            print 'Create node vector'
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

            UDN[mo,:] = uDN_node/ind_cnt
            VDN[mo,:] = vDN_node/ind_cnt
            UVDN[mo,:] = uvDN_node/ind_cnt

        tracershape = np.shape(UDN)

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
          w_nc_fid.description = u'HEDY' 
          w_nc_fid.history     = 'Created ' + time.ctime(time.time())

          nod2d    = w_nc_fid.createDimension('nod2d', mesh.n2d)               # Create dimension: number of 3d nodes
          time	   = w_nc_fid.createDimension('time', tracershape[0]) 

          w_nc_var = w_nc_fid.createVariable('HEDY', 'f4',('time','nod2d'))           # 'DIN' is name of saved variable/ 'f8' sets presicion to 64-bit floating point
          w_nc_var.setncatts({'long_name': u'HEDY',                          'units': u'mmol/m2/day'})
          w_nc_fid.variables['HEDY'][:] = UVDN   

          w_nc_var = w_nc_fid.createVariable('UEDY', 'f4',('time','nod2d'))           # 'DIN' is name of saved variable/ 'f8' sets presicion to 64-bit floating point
          w_nc_var.setncatts({'long_name': u'UEDY',                          'units': u'mmol/m2/day'})
          w_nc_fid.variables['UEDY'][:] = UDN  

          w_nc_var = w_nc_fid.createVariable('VEDY', 'f4',('time','nod2d'))           # 'DIN' is name of saved variable/ 'f8' sets presicion to 64-bit floating point
          w_nc_var.setncatts({'long_name': u'VEDY',                          'units': u'mmol/m2/day'})
          w_nc_fid.variables['VEDY'][:] = VDN  

          w_nc_fid.close()                                                     # close the new file                

          cwd = os.getcwd()
          print "New netcdf file (",netcdf_name,") has been created."
          print "Location: "+output_directory
        else:
          print 'You have specified not to save your field in netcdf file'
    

edyhv2(sys.argv[1])