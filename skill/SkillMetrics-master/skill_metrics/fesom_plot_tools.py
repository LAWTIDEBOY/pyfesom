import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
def ftriplot(mesh, data2, contours, cmap=[], oce='global', do_cbar=True, mlabels=[0,0,0,0], plabels=[0,0,0,0], extend='both'):
	if (cmap==[]):
		cmap=plt.cm.jet
	if (oce=='global'):
		data2=np.copy(data2)
		elem2=mesh.elem[mesh.no_cyclic_elem,:]
		d=data2[elem2].mean(axis=1)
		k = [i for (i, val) in enumerate(d) if not np.isnan(val)]
		elem2=elem2[k,:]
		print 'ftriplot, number of dummy points:', len(d)-len(k)	
		map = Basemap(projection='robin',lon_0=0)
		x, y = map(mesh.x2, mesh.y2)
		map.drawmapboundary(fill_color='0.9')
		map.drawcoastlines()
		map.drawparallels(np.arange(-90,90,30),labels=plabels) #[1,0,0,0]
		map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=mlabels) #[0,0,0,1]
		#data2[data2>900]=np.nan
		eps=(contours.max()-contours.min())/50.
		data2[data2<=contours.min()]=contours.min()+eps
		data2[data2>=contours.max()]=contours.max()-eps
		im=plt.tricontourf(x, y, elem2, data2, levels=contours, cmap=cmap, extend=extend)
		if do_cbar:
			cbar=map.colorbar(im,"bottom", size="5%", pad="2%")
	elif (oce=='np'):
		data2=np.copy(data2)
		elem2=mesh.elem#[mesh.no_cyclic_elem,:]
		d=data2[elem2].mean(axis=1)
		k = [i for (i, val) in enumerate(d) if not np.isnan(val)]
		elem2=elem2[k,:]
		print 'ftriplot, number of dummy points:', len(d)-len(k)		
	        map = Basemap(projection='nplaea',boundinglat=45,lon_0=0,resolution='l')
		x, y = map(mesh.x2, mesh.y2)
		map.drawcoastlines()
		map.drawparallels(np.arange(-80.,81.,20.), labels=plabels)
		map.drawmeridians(np.arange(-180.,181.,20.),labels=mlabels) #[0,1,0,0]
		map.drawmapboundary(fill_color='0.9')
		map.fillcontinents(color='.7',lake_color='.7')
		#data2[data2>900]=np.nan
		eps=(contours.max()-contours.min())/100.
		data2[data2<=contours.min()]=contours.min()+eps
		data2[data2>=contours.max()]=contours.max()-eps
		im=plt.tricontourf(x, y, elem2, data2, levels=contours, cmap=cmap, extend=extend)
		if do_cbar:
			cbar=map.colorbar(im,"bottom", size="5%", pad="2%")
	elif (oce=='sp'):
		data2=np.copy(data2)
		elem2=mesh.elem#[mesh.no_cyclic_elem,:]
		d=data2[elem2].mean(axis=1)
		k = [i for (i, val) in enumerate(d) if not np.isnan(val)]
		elem2=elem2[k,:]
		print 'ftriplot, number of dummy points:', len(d)-len(k)		
		map = Basemap(projection='splaea',boundinglat=-20,lon_0=180,resolution='l')
		x, y = map(mesh.x2, mesh.y2)
		map.drawcoastlines()
		map.drawparallels(np.arange(-80.,81.,20.), labels=plabels)
		map.drawmeridians(np.arange(-180.,181.,20.),labels=mlabels)
		map.drawmapboundary(fill_color='0.9')
		map.fillcontinents(color='.7',lake_color='.7')
		#data2[data2>900]=np.nan
		eps=(contours.max()-contours.min())/100.
		data2[data2<=contours.min()]=contours.min()+eps
		data2[data2>=contours.max()]=contours.max()-eps
		im=plt.tricontourf(x, y, elem2, data2, levels=contours, cmap=cmap, extend=extend)
		if do_cbar:
			cbar=map.colorbar(im,"bottom", size="5%", pad="2%")
			cbar.set_ticks([round(i,1) for i in np.linspace(cont[0], cont[1], 5)])
	return(im, cbar if (do_cbar) else False)

def wplot_xy(xx,yy,zz,contours, cmap=[]):
	import numpy as np
	import matplotlib.pyplot as plt
	from mpl_toolkits.basemap import Basemap
	from matplotlib.colors import LinearSegmentedColormap

	if (cmap==[]):
		cmap=plt.cm.jet
	eps=(contours.max()-contours.min())/100.
	zz[zz<=contours.min()]=contours.min()+eps
	zz[zz>=contours.max()]=contours.max()-eps

	map = Basemap(projection='robin',lon_0=0)
	xxx, yyy = map(xx, yy)

	map.drawmapboundary(fill_color='0.9')
	map.drawcoastlines()
	map.fillcontinents(color='.7',lake_color='.7')
	map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
	map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=[0,0,0,1])

	im=plt.contourf(xxx, yyy, zz, levels=contours, cmap=cmap, extend='both')
	cbar=map.colorbar(im,"bottom", size="5%", pad="2%")
	return(im, cbar)

def wplot_yz(y,z,v,contours, cmap=[]):
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.colors import LinearSegmentedColormap

	if (cmap==[]):
		cmap=plt.cm.jet

	im=plt.contourf(y, z, v, levels=contours, cmap=cmap, extend='both')
	cbar=plt.colorbar(orientation='horizontal')
	plt.grid()
	return(im, cbar)

def movingaverage(interval, window_size):
	import numpy as np
	window= np.ones(int(window_size))/float(window_size)	
	ret=list(interval)
	for i in range (window_size):
		ret=ret+[ret[-1]]
	ret=np.convolve(np.array(ret), window, 'valid')
	return ret

