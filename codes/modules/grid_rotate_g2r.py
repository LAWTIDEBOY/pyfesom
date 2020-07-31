def grid_rotate_g2r(al,be,ga,lon,lat):
    import numpy as np
    import math as mt

    rad=mt.pi/180
    al=al*rad
    be=be*rad
    ga=ga*rad

    rotate_matrix=np.zeros(shape=(3,3))
    rotate_matrix[0,0]=np.cos(ga)*np.cos(al)-np.sin(ga)*np.cos(be)*np.sin(al)
    rotate_matrix[0,1]=np.cos(ga)*np.sin(al)+np.sin(ga)*np.cos(be)*np.cos(al)
    rotate_matrix[0,2]=np.sin(ga)*np.sin(be)
    rotate_matrix[1,0]=-np.sin(ga)*np.cos(al)-np.cos(ga)*np.cos(be)*np.sin(al)
    rotate_matrix[1,1]=-np.sin(ga)*np.sin(al)+np.cos(ga)*np.cos(be)*np.cos(al)
    rotate_matrix[1,2]=np.cos(ga)*np.sin(be)
    rotate_matrix[2,0]=np.sin(be)*np.sin(al)
    rotate_matrix[2,1]=-np.sin(be)*np.cos(al)
    rotate_matrix[2,2]=np.cos(be)

    lat=lat*rad
    lon=lon*rad

    xr=np.cos(lat)*np.cos(lon)
    yr=np.cos(lat)*np.sin(lon)
    zr=np.sin(lat)	#Geographical Cartesian coordinates:

    xg=rotate_matrix[0,0]*xr + rotate_matrix[0,1]*yr + rotate_matrix[0,2]*zr
    yg=rotate_matrix[1,0]*xr + rotate_matrix[1,1]*yr + rotate_matrix[1,2]*zr
    zg=rotate_matrix[2,0]*xr + rotate_matrix[2,1]*yr + rotate_matrix[2,2]*zr		

    rlat=[mt.asin(val) for (i, val) in enumerate(zg)]
    rlon=[mt.atan2(yg[i],xg[i]) for i in range(len(xg))]
    a = [i for (i, val) in enumerate((abs(xg)+abs(yg))) if val ==0]
    #.astype(float32)
    if a: rlon[a]=0
    rlat = [rlat[i]/rad for i in range(len(rlat))]
    rlon = [rlon[i]/rad for i in range(len(rlon))]
    rlat=np.array(rlat)
    rlon=np.array(rlon)

    return (rlon,rlat)
