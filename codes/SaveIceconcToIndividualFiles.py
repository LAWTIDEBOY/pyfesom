import sys 
sys.path.append('../../modules')

import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from load_mesh_data_new import load_mesh
import numpy as np
from netCDF4 import Dataset
from PiecewiseNorm import PiecewiseNorm
import os.path
import time

# ==============================================================================
# Running this file loads tracers from old FESOM-REcoM2 output file (oce.mean.nc)
# and saves each tracer in an individual file
# 
#  Input:
#  - mesh_id: Name of mesh, will be added to the netcdf name
#  - meshpath: Speciefies where the target mesh is stored
#  - save_netcdf: If true, netcdf will be created
#  - delete_old_netcdf: If a netcdf file with the same name exists, a new cannot 
#    be made. If set to true, an old netcdf with the same name will be deleted
#
#  Output:
#  - netcdf file for each tracer in the old file
#  
#  During running, keep an eye on the output in the terminal, to see if it 
#  makes sense. 
#
# ==============================================================================
# Loading mesh for run

mesh_id    = 'meshArc4.5'
meshpath   = '/scratch/usr/hbkvsk12/hlrn3_work2/mesh/'+mesh_id+'/'            # Defining path where mesh is stored
mesh       = load_mesh(meshpath)                                    # Loading mesh, stores it in mesh.****  
mesh.zlevs = -mesh.zlevs                                            # Depth is made negative

tracername = 'area'
year	= 2015
runid	= 'Arc12'

# ==============================================================================
# Settings for netcdf file

save_netcdf       = True                                            # Saves the interpolated field in netcdf file
delete_old_netcdf = True                                            # If a netcdf file with the same name exists it will be deleted
saving_directory  = '/scratch/usr/hbkvsk12/hlrn3_work2/results/Arc12/'                                    # Where the netcdf is saved
netcdf_name       = runid +'.'+str(year)+'.monthly.ice.mean.nc'                        # Name of netcdf with interpolated DIN field
plot_netcdf       = False                                            # Reads DIN from the created netcdf file, else it plots the interpolated field (should be the same)

# ==============================================================================
# Loading data

ncfile	= '/scratch/usr/hbkvsk12/hlrn3_work2/results/'+runid+'/Oldfiles/'+runid+'.'+str(year)+'.ice.mean.nc'
f	= Dataset(ncfile,'r')
tracer	= f.variables['area'][:]
t1 = np.mean(tracer[0:15,:],axis=0)
t2 = np.mean(tracer[15:29,:],axis=0)
t3 = np.mean(tracer[29:45,:],axis=0)
t4 = np.mean(tracer[45:60,:],axis=0)
t5 = np.mean(tracer[60:75,:],axis=0)
t6 = np.mean(tracer[75:90,:],axis=0)
t7 = np.mean(tracer[90:106,:],axis=0)
t8 = np.mean(tracer[106:121,:],axis=0)
t9 = np.mean(tracer[121:136,:],axis=0)
t10 = np.mean(tracer[136:152,:],axis=0)
t11 = np.mean(tracer[152:167,:],axis=0)
t12 = np.mean(tracer[166:,:],axis=0)
tracer = [t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12]

tracershape = np.shape(tracer)

# ==============================================================================
# Testing if a netcdf file with the same name exists, if yes, it must be removed
# to save a new one.

if os.path.isfile(saving_directory+netcdf_name) and delete_old_netcdf:
  os.remove(saving_directory+netcdf_name)
  print "The netcdf file "+netcdf_name+" has been deleted to make room for your file of the same name."
elif os.path.isfile(netcdf_name):
  statement = "The netcdf file "+netcdf_name+" already exists! It must be removed for a new one to be created. This can be done by changing your settings."
  sys.exit(statement)
  
if not os.path.isdir(saving_directory):
  os.makedirs(saving_directory)
  print 'Directory '+saving_directory+' has been created'
  
# ==============================================================================
# Creating netcdf file
if save_netcdf:    
  w_nc_fid = Dataset(saving_directory+netcdf_name, 'w', format='NETCDF4_CLASSIC')      # Create and open new netcdf file to write to
  w_nc_fid.description = "Sea ice conc" 
  w_nc_fid.history     = 'Created ' + time.ctime(time.time())

  nod3d    = w_nc_fid.createDimension('nod2d', mesh.n2d)               # Create dimension: number of 3d nodes
  time	   = w_nc_fid.createDimension('time', tracershape[0]) 

  w_nc_var = w_nc_fid.createVariable(tracername, 'f8',('time','nod2d'))           # 'DIN' is name of saved variable
                                                                       # 'f8' sets presicion to 64-bit floating point
  w_nc_var.setncatts({'long_name': u"Sea ice conc",\
                      'units': u"0 to 1"})
  w_nc_fid.variables[tracername][:] = tracer   
  w_nc_fid.close()                                                     # close the new file                
  
  cwd = os.getcwd()
  print "New netcdf file (",netcdf_name,") has been created."
  print "Location: "+saving_directory
else:
  print 'You have specified not to save your field in netcdf file'  
  
# ==============================================================================
# Plot the surface of interpolated DIN field

print '***' 
if plot_netcdf:
  ncfile     = saving_directory+netcdf_name
  f          = Dataset(ncfile, 'r')
  data1      = f.variables[tracername][0,0:len(mesh.x2)]          
  print "Plotting alk from netcdf file"
else:  
  data1 = din_int[0:len(mesh.x2)]
  print "Plotting interpolated phy N field"
  
contours = (0.,1.,.01)
contours = np.arange(contours[0], contours[1]+contours[2], contours[2])

fig = plt.figure(num=1, figsize=(12, 8), facecolor='w', edgecolor='k')
ax=plt.subplot(111)
elem2=mesh.elem[mesh.no_cyclic_elem,:]
d=data1[elem2].mean(axis=1)
k = [i for (i, val) in enumerate(d) if not np.isnan(val)]
elem2=elem2[k,:]    
bmap = Basemap(projection='robin',lon_0=0) 
x, y = bmap(mesh.x2, mesh.y2)
bmap.drawmapboundary(fill_color='0.9')
bmap.drawcoastlines()
mlabels=[False,False,False,False]
plabels=[True,True,True,True]		    
bmap.drawparallels(np.arange(-90,90,30),labels=plabels) 
bmap.drawmeridians(np.arange(bmap.lonmin,bmap.lonmax+30,60),labels=mlabels)
im=plt.tricontourf(x, y, elem2, data1, levels=contours, cmap='rainbow', norm=PiecewiseNorm(contours))
label = 'Surface PhyN'
plt.title(label)
cbar=bmap.colorbar(im,"bottom", size="5%", pad="2%")
cbar.set_label(r'[mmol m$^{-3}$]')

#plt.savefig(saving_directory+'ChlNano.png', dpi = 200, bbox_inches='tight')
plt.show()
