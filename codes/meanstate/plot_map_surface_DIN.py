{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 77,
   "metadata": {},
   "outputs": [],
   "source": [
    "import os\n",
    "import sys\n",
    "sys.path.append('modules')\n",
    "\n",
    "import matplotlib.pyplot as plt\n",
    "from mpl_toolkits.basemap import Basemap\n",
    "from pylab import *\n",
    "from load_mesh_data_new import *\n",
    "import numpy as np\n",
    "from PiecewiseNorm import PiecewiseNorm\n",
    "from netCDF4 import Dataset\n",
    "import colormaps as cmaps\n",
    "from matplotlib.colors import ListedColormap\n",
    "import earthpy as et "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 78,
   "metadata": {},
   "outputs": [],
   "source": [
    "my_path=os.path.join(et.io.HOME,\"Documents\",\"projects\",\"awi-models\")\n",
    "#os.path.exists(my_path)\n",
    "os.chdir(my_path)\n",
    "#os.getcwd()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 79,
   "metadata": {},
   "outputs": [],
   "source": [
    "viridis = ListedColormap(cmaps.viridis.colors)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 80,
   "metadata": {},
   "outputs": [],
   "source": [
    "global_plot = False\n",
    "arctic_plot = True"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 81,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "'data/meshArc4.5/'"
      ]
     },
     "execution_count": 81,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "runid='ArcC3'\n",
    "first_year = 1980\n",
    "last_year  = 1980\n",
    "years      = np.arange(first_year,last_year+1,1)\n",
    "months =np.linspace(0,11,12).astype(int)\n",
    "\n",
    "resultpath = 'data/'+runid+'/'\n",
    "savepath    = 'figures/'+runid+'/'\n",
    "meshpath    = 'data/meshArc4.5/'"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 90,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "The *usepickle = True* and the pickle file (*pickle_mesh*) exists.\n",
      " We load the mesh from it.\n"
     ]
    },
    {
     "ename": "TypeError",
     "evalue": "a bytes-like object is required, not 'str'",
     "output_type": "error",
     "traceback": [
      "\u001b[0;31m---------------------------------------------------------------------------\u001b[0m",
      "\u001b[0;31mTypeError\u001b[0m                                 Traceback (most recent call last)",
      "\u001b[0;32m<ipython-input-90-dc1d635fc790>\u001b[0m in \u001b[0;36m<module>\u001b[0;34m\u001b[0m\n\u001b[1;32m      1\u001b[0m \u001b[0;31m#mesh = fesom_mesh(meshpath, get3d = False)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m----> 2\u001b[0;31m \u001b[0mmesh\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mload_mesh\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mmeshpath\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m",
      "\u001b[0;32m~/Documents/projects/awi-models/codes/modules/load_mesh_data_new.py\u001b[0m in \u001b[0;36mload_mesh\u001b[0;34m(path, abg, get3d, usepickle)\u001b[0m\n\u001b[1;32m     30\u001b[0m \u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     31\u001b[0m         \u001b[0mifile\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mopen\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mpickle_file\u001b[0m\u001b[0;34m,\u001b[0m \u001b[0;34m'r'\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0;32m---> 32\u001b[0;31m         \u001b[0mmesh\u001b[0m \u001b[0;34m=\u001b[0m \u001b[0mpickle\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mload\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0mifile\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[0m\u001b[1;32m     33\u001b[0m         \u001b[0mifile\u001b[0m\u001b[0;34m.\u001b[0m\u001b[0mclose\u001b[0m\u001b[0;34m(\u001b[0m\u001b[0;34m)\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n\u001b[1;32m     34\u001b[0m         \u001b[0;32mreturn\u001b[0m \u001b[0mmesh\u001b[0m\u001b[0;34m\u001b[0m\u001b[0;34m\u001b[0m\u001b[0m\n",
      "\u001b[0;31mTypeError\u001b[0m: a bytes-like object is required, not 'str'"
     ]
    }
   ],
   "source": [
    "#mesh = fesom_mesh(meshpath, get3d = False)\n",
    "mesh = load_mesh(meshpath)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {},
   "outputs": [],
   "source": [
    "# load the given biological tracer #\n",
    "var_id = 'tr01' #DIN"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 88,
   "metadata": {},
   "outputs": [
    {
     "ename": "SyntaxError",
     "evalue": "Missing parentheses in call to 'print'. Did you mean print(years[ind])? (<ipython-input-88-42d5abf8b87a>, line 3)",
     "output_type": "error",
     "traceback": [
      "\u001b[0;36m  File \u001b[0;32m\"<ipython-input-88-42d5abf8b87a>\"\u001b[0;36m, line \u001b[0;32m3\u001b[0m\n\u001b[0;31m    print years[ind]\u001b[0m\n\u001b[0m              ^\u001b[0m\n\u001b[0;31mSyntaxError\u001b[0m\u001b[0;31m:\u001b[0m Missing parentheses in call to 'print'. Did you mean print(years[ind])?\n"
     ]
    }
   ],
   "source": [
    "DIN2D = np.zeros(len(mesh.x2))\n",
    "for ind in range(0,len(years)):\n",
    "    print years[ind]\n",
    "    ncfile = resultpath+runid+'.'+str(years[ind])+'.oce.mean.nc'\n",
    "    f      = Dataset(ncfile, 'r')    \n",
    "    din    = f.variables[var_id][:].mean(axis=0)\n",
    "    DIN2D  = DIN2D + din[0:mesh.n2d]   \n",
    "data2 = DIN2D/len(years)\n",
    "\n",
    "print 'Number of nans in tracer: ',var_id,np.count_nonzero(np.isnan(data2))\n",
    "print 'Number of inf in tracer: ',var_id,np.count_nonzero(np.isinf(data2))\n",
    "print 'Mean of surface: ',np.mean(data2[0:len(mesh.x2)])\n",
    "print 'Max and min: ',np.max(data2),np.min(data2)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 89,
   "metadata": {},
   "outputs": [
    {
     "ename": "SyntaxError",
     "evalue": "invalid syntax (<ipython-input-89-4d59192f923f>, line 1)",
     "output_type": "error",
     "traceback": [
      "\u001b[0;36m  File \u001b[0;32m\"<ipython-input-89-4d59192f923f>\"\u001b[0;36m, line \u001b[0;32m1\u001b[0m\n\u001b[0;31m    print 'Number of nans in tracer: ',var_id,np.count_nonzero(np.isnan(data2))\u001b[0m\n\u001b[0m                                     ^\u001b[0m\n\u001b[0;31mSyntaxError\u001b[0m\u001b[0;31m:\u001b[0m invalid syntax\n"
     ]
    }
   ],
   "source": [
    "print 'Number of nans in tracer: ',var_id,np.count_nonzero(np.isnan(data2))\n",
    "print 'Number of inf in tracer: ',var_id,np.count_nonzero(np.isinf(data2))\n",
    "print 'Mean of surface: ',np.mean(data2[0:len(mesh.x2)])\n",
    "print 'Max and min: ',np.max(data2),np.min(data2)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "contours = [0.0, 8., 1.]\n",
    "    contours = np.arange(contours[0], contours[1]+contours[2], contours[2])\n",
    "    fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')\n",
    "    elem2=mesh.elem#[mesh.no_cyclic_elem,:]\n",
    "    d=data2[elem2].mean(axis=1)\n",
    "    k = [i for (i, val) in enumerate(d) if not np.isnan(val)]\n",
    "    elem2=elem2[k,:]\n",
    "    print 'ftriplot, number of dummy points:', len(d)-len(k)\t\t\n",
    "    map = Basemap(projection='nplaea',boundinglat=65,lon_0=0,resolution='l')\n",
    "    x, y = map(mesh.x2, mesh.y2)\n",
    "    map.drawcoastlines()\n",
    "    plabels=[False,False,False,False]\n",
    "    mlabels=[True,True,True,True]    \n",
    "    map.drawparallels(np.arange(-80.,81.,20.), labels=plabels)\n",
    "    map.drawmeridians(np.arange(-180.,181.,20.),labels=mlabels) #[0,1,0,0]\n",
    "    map.drawmapboundary(fill_color='0.9')\n",
    "    map.fillcontinents(color='.7',lake_color='.7')\n",
    "    eps=(contours.max()-contours.min())/100.\n",
    "    data2[data2<=contours.min()]=contours.min()+eps\n",
    "    data2[data2>=contours.max()]=contours.max()-eps\n",
    "    im=plt.tricontourf(x, y, elem2, data2, levels=contours, cmap=viridis, norm=PiecewiseNorm(contours))\n",
    "    label = runid+': Mean surface DIN conc. ('+str(years[0])+' to '+str(years[len(years)-1])+')'\n",
    "    plt.title(label,y=1.05)\n",
    "    cbar=map.colorbar(im,\"bottom\", size=\"5%\", pad=\"5%\")\n",
    "    cbar.set_label(r'DIN conc. [mmol N m$^{-3}$]')   \n",
    "    plt.savefig(savepath+'DINspatialArc'+str(years[ind])+'.png', dpi = dpicnt, bbox_inches='tight') \n",
    "\n",
    "plt.show()   "
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
