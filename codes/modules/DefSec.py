def DefSec(num_section,mesh,pos):
  import numpy as np
  import math as mt
  print 'Searching for section'	

  reso = 10

  xsur = mesh.x2[0:mesh.n2d]
  ysur = mesh.y2[0:mesh.n2d]

  nod32 = mesh.n32.T-1
  elem=mesh.elem[mesh.no_cyclic_elem,:]
  elem=np.matrix.transpose(elem)

  k_all, m_all = np.zeros(num_section), np.zeros(num_section)
  n1_all, n2_all, intcoef_all, vol_el_all = {}, {}, {}, {}
  nvec_all = np.zeros(shape=(3,num_section))
  n_cnt, tricount_cnt = np.zeros(num_section), np.zeros(num_section)

  for nsec in range(0,num_section):  
    xtmp=[]
    ytmp=[]
  
    x1, y1, x2, y2 = pos[0,nsec], pos[1,nsec], pos[2,nsec], pos[3,nsec]
    dirx=(x2-x1)/np.sqrt((x2-x1)**2+(y2-y1)**2)
    diry=(y2-y1)/np.sqrt((x2-x1)**2+(y2-y1)**2)
    tri=[]
    section = np.empty((0,2), int)
    intcoef=[]
    dst=[]
    ind2d=mesh.nodes
    bnod = [i for i, x in enumerate(ind2d) if x == 1]
    bnod = np.array(bnod)
    bn_x, bn_y = xsur[bnod], ysur[bnod]
    a, b = (y1-y2)/(x1-x2), -1 
    c = y1-a*x1 
    d = (a*bn_x+b*bn_y+c)/np.sqrt(a**2+b**2) 

    a1 = np.sqrt(np.cos(y1/180*mt.pi)**2*(bn_x-x1)**2+(bn_y-y1)**2)*mt.pi/180.*6367.5
    a1 = ((d < 0) & (a1 < 4*reso)).nonzero()[0]
    a1 = np.array(a1)
    n=np.array([bnod[a1]])
    maxval = max(d[a1])
    a2 = np.where(d[a1] == d[a1].max())
    node1=n[0,a2]
  
    t=ind2d[elem[0,:]]+ind2d[elem[1,:]]+ind2d[elem[2,:]]
    ee = ((np.sum(1*(elem==node1),axis=0)==1) & (t==2)).nonzero()[0]

    nn = len(ee)
    flag = 0
    for i in range(0,nn):
      elnodes=elem[:,ee[i]]
      a2 = (a*xsur[elnodes]+b*ysur[elnodes]+c)/np.sqrt(a**2+b**2)
      a2 = ((a2 > 0) & (ind2d[elnodes] ==1)).nonzero()[0]
      if a2.size:
        flag=1
        break

    if flag==0:
      a1 = np.sqrt(np.cos(y1/180*mt.pi)**2*(bn_x-x1)**2+(bn_y-y1)**2)*mt.pi/180*6367.5
      a1 = ((d > 0) & (a1 < 10*reso)).nonzero()[0]
      n  = bnod[a1]
      maxval = min(d[a1])
      a2 = np.where(d[a1] == d[a1].min())
      node1=n[a2]
      t=ind2d[elem[0,:]]+ind2d[elem[1,:]]+ind2d[elem[2,:]]
      ee = ((np.sum(1*(elem==node1),axis=0)==1) & (t==2)).nonzero()[0]
      nn = len(ee)

      flag=0
      for i in range(0,nn):
        elnodes=elem[:,ee[i]]
        a2 = (a*xsur[elnodes]+b*ysur[elnodes]+c)/np.sqrt(a**2+b**2)
        a2 = ((a2 > 0) & (ind2d[elnodes] ==1)).nonzero()[0]
        if a2.size:
          flag=1
          break

    if flag==0:
      print 'Cannot define the section line.'
      sys.exit()

    node2   = elnodes[a2]
    row     = np.array([node1,node2],dtype=np.int64).T
    section = np.append(section,row,axis=0)#np.vstack((section, row))
    xp, yp  = xsur[node1], ysur[node1]
    betax, betay = xsur[node2]-xp, ysur[node2]-yp
    d       = -diry*betax+dirx*betay
    d1      = (diry*(xp-x1)-dirx*(yp-y1))/d
    c1      = (betay*(xp-x1)-betax*(yp-y1))/d
    intcoef = np.append(intcoef, np.array([[d1]]))
    dst     = np.append(dst, np.array([[c1]]))
    xtmp    = np.append(xtmp, np.array([[xp+betax*intcoef[0]]]))
    ytmp    = np.append(ytmp, np.array([[yp+betay*intcoef[0]]]))

    # Now the regular cycle to find all other points:
    tricount = 0
    success  = 0
    elem0    = elem

    while success==0:
      node1=section[tricount,0]
      node2=section[tricount,1]
      # Find third node that completes triangles. There should be 1 or 2 such nodes
      a1    = ((elem0[0,:]==node1) | (elem0[1,:]==node1) | (elem0[2,:]==node1)).nonzero()[0]
      elem1 = elem0[:,a1]
      a2    = ((elem1[0,:]==node2) | (elem1[1,:]==node2) | (elem1[2,:]==node2)).nonzero()[0]
      if len(a2)==1 or (nsec==7 and xtmp[-1]<3) or (nsec==6 and xtmp[-1]<-5):
        if tricount>1:
          success=1
        else:
          tricount += 1  
          tri       = np.append(tri, np.array([a1[a2[0]]]))
      else:
        tri1=a1[a2[0]]
        tri2=a1[a2[1]]
        tricount += 1
        if tri1==tri[tricount-2]:
          tri = np.hstack((tri,tri2))
        else:
          tri = np.hstack((tri,tri1)) 
      
      if success==0:
        # tri(tricount) contains the number of the proper triangle. 
        elnodes=elem0[:,tri[tricount-1]]
        a2=((elnodes!=node1) & (elnodes!=node2))
        node3=elnodes[a2]
        # First possibility:
        xp, yp = xsur[node1], ysur[node1]
        xf, yf = xsur[node3], ysur[node3]
        betax, betay = xf-xp, yf-yp
        d     = -diry*betax+dirx*betay
        if d==0:
          d=0.00000001
        d1 = (diry*(xp-x1)-dirx*(yp-y1))/d
        c1 = (betay*(xp-x1)-betax*(yp-y1))/d
        if d1>=0 and d1<=1:
          row     = np.array([node1,node3],dtype=np.int64).T
          section = np.append(section,[row],axis=0)
          intcoef = np.append(intcoef, np.array([[d1]]))
          dst     = np.append(dst, np.array([[c1]])) 
        else:
          # Second possibility:
          xp, yp = xsur[node2], ysur[node2]
          betax, betay = xf-xp, yf-yp
          d     = -diry*betax+dirx*betay
          if d==0: 
            d=0.00000001
          d2=(diry*(xp-x1)-dirx*(yp-y1))/d
          c2=(betay*(xp-x1)-betax*(yp-y1))/d
          # intersection requires d between 0 and 1 
          row     = np.array([node2,node3],dtype=np.int64).T
          section = np.append(section,[row],axis=0)
          intcoef = np.append(intcoef, np.array([[d2]]))
          dst     = np.append(dst, np.array([[c2]])) 
        xtmp    = np.append(xtmp, np.array([[xp+betax*intcoef[tricount]]]))
        ytmp    = np.append(ytmp, np.array([[yp+betay*intcoef[tricount]]]))  

    n1 = nod32[:,section[:,0]]
    n2 = nod32[:,section[:,1]]   
    levels=np.zeros(tricount+1)
    for nn in range(0,tricount+1):
      a1 = np.nonzero((n1[:,nn] == -1000)) # find index of non-existing values
      if not a1:
        a1 = np.shape(a1)+1
      a2 = np.nonzero((n2[:,nn] == -1000))   
      if not a2:
        a2 = np.shape(a2)+1
        
      levels[nn] = a1[0][0]-1
      
      n2[a2[0][0]:len(n2),nn] = n2[a2[0][0]-1,nn]
      n1[a1[0][0]:len(n1),nn] = n1[a1[0][0]-1,nn]
     
      if a1[0][0]>a2[0][0]:
        levels[nn] = a2[0][0]
      if a2[0][0]>a1[0][0]:
        levels[nn] = a1[0][0]    
  
    # tricount+1 is the length of array sect
    n=np.max(levels)+1
    zs=np.empty(shape=(n, tricount+1))
    zs[:] = np.NAN
    xs=np.zeros(shape=(n, tricount+1))
    ss=np.zeros(shape=(n, tricount+1))
    ts=np.zeros(shape=(n, tricount+1))
    us=np.zeros(shape=(n, tricount+1))
    vs=np.zeros(shape=(n, tricount+1))
    
    for nn in range(0,tricount+1):
      nk = levels[nn]+1
      zs[0:nk,nn] = mesh.z3[n1[0:nk,nn]]*(1-intcoef[nn])+mesh.z3[n2[0:nk,nn]]*intcoef[nn] 
      xs[:,nn]    = dst[nn]
     
    # area
    dz1 = zs[0:len(zs)-1,:]-zs[1:len(zs),:] 
    dz  = (dz1[:,0:len(zs[0,:])-1]+dz1[:,1:len(zs[0,:])])/2.0  
    
    lon = xtmp
    lat = ytmp
   
    earth  = 6367.5e3  #fesom parameter
    deg2rad= mt.pi/180
    nx     = len(lon)
    lt = deg2rad*lat
    ln = deg2rad*lon
    alfa = np.arccos( np.cos(lt[0:nx-1])*np.cos(lt[1:nx])*np.cos(ln[0:nx-1]-ln[1:nx])+np.sin(lt[0:nx-1])*np.sin(lt[1:nx]))
    dx = earth*abs(alfa)
    vol_el = np.tile(dx,[n-1,1])*dz
    
    # normal direction
    vec1 = np.array([0, 0, 1])
    vec2 = np.array([(x2-x1)*np.cos((y2+y1)/2*mt.pi/180), y2-y1, 0])
    vec2 = vec2/np.sqrt(vec2.dot(vec2.T))
    nv_se = np.zeros(3)
    nv_se[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1]
    nv_se[1]=vec1[2]*vec2[0]-vec1[0]*vec2[2]
    nv_se[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0]
    nvec = -nv_se/np.sqrt(nv_se.dot(nv_se.T))
    
    # Interpolation coefficient
    intcoef = np.tile(intcoef,[n,1])
    n1=n1[0:n,:]
    n2=n2[0:n,:]
    
    # longitude
    lon = (xtmp[0:len(zs[0,:])-1]+xtmp[1:len(zs[0,:])])/2.0
    lon = np.tile(lon,[n-1,1])
    
    # save
    k_all[nsec]=n
    m_all[nsec]=tricount+1
    nvec_all[:,nsec] = nvec
    for ind1 in range(0,int(n)): # A dictionary with all values are created, this must be read and added to matrix before it can be used as index
      for ind2 in range(0,int(tricount)+1):
        n1_all[ind1,ind2,nsec]      = n1[ind1,ind2]
        n2_all[ind1,ind2,nsec]      = n2[ind1,ind2]
        intcoef_all[ind1,ind2,nsec] = intcoef[ind1,ind2]
    for ind1 in range(0,int(n)-1): # A dictionary with all values are created, this must be read and added to matrix before it can be used as index
      for ind2 in range(0,int(tricount)):
        vol_el_all[ind1,ind2,nsec]  = vol_el[ind1,ind2]
        
    n_cnt[nsec]        = n
    tricount_cnt[nsec] = tricount
    print 'Found nodes for section', nsec

  # Recreating matrixes saved in tuples  
  n = int(np.max(n_cnt))
  tricount = int(np.max(tricount_cnt))

  n1_all1      = np.zeros(shape=(n,tricount+1,num_section))
  n2_all1      = np.zeros(shape=(n,tricount+1,num_section))
  intcoef_all1 = np.zeros(shape=(n,tricount+1,num_section))
  vol_el_all1  = np.zeros(shape=(n-1,tricount,num_section))
 
  for nsec in range(0,num_section):
    for ind1 in range(0,int(n_cnt[nsec])):
      for ind2 in range(0,int(tricount_cnt[nsec]+1)):
        n1_all1[ind1,ind2,nsec]      = n1_all[ind1,ind2,nsec]
        n2_all1[ind1,ind2,nsec]      = n2_all[ind1,ind2,nsec]
        intcoef_all1[ind1,ind2,nsec] = intcoef_all[ind1,ind2,nsec] 
    for ind1 in range(0,int(n_cnt[nsec])-1):
      for ind2 in range(0,int(tricount_cnt[nsec])):
        vol_el_all1[ind1,ind2,nsec]  = vol_el_all[ind1,ind2,nsec]   

  n1_all, n2_all, intcoef_all, vol_el_all = n1_all1.astype(int), n2_all1.astype(int), intcoef_all1, vol_el_all1   

  return (n1_all,n2_all,intcoef_all,vol_el_all,k_all,m_all,nvec_all)


