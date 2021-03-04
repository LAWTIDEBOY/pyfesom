def vec_rotate_r2g(lon, lat,tlon, tlat):
    import numpy as np
    import math as mt
    from rotation import grid_rotate_g2r

    al, be, ga = 50, 15, -90

    (rlon,rlat) = grid_rotate_g2r(al,be,ga,lon,lat)
 
    result = []
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

    rotate_matrix=np.linalg.pinv(rotate_matrix)

    rlat=rlat*rad
    rlon=rlon*rad	#Rotated Cartesian coordinates:
    lat = lat*rad
    lon=lon*rad

    txg=-tlat*np.sin(rlat)*np.cos(rlon)-tlon*np.sin(rlon)
    tyg=-tlat*np.sin(rlat)*np.sin(rlon)+tlon*np.cos(rlon)
    tzg= tlat*np.cos(rlat)	#Geographical Cartesian coordinates:

    txr=rotate_matrix[0,0]*txg + rotate_matrix[0,1]*tyg + rotate_matrix[0,2]*tzg
    tyr=rotate_matrix[1,0]*txg + rotate_matrix[1,1]*tyg + rotate_matrix[1,2]*tzg
    tzr=rotate_matrix[2,0]*txg + rotate_matrix[2,1]*tyg + rotate_matrix[2,2]*tzg		#Geographical coordinates:
    v=-np.sin(lat)*np.cos(lon)*txr - np.sin(lat)*np.sin(lon)*tyr + np.cos(lat)*tzr
    u=-np.sin(lon)*txr + np.cos(lon)*tyr

    return (u,v)
