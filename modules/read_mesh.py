class fesom_mesh:
  """existing instances are: path, n2d, e2d, nlev, zlevs, x2, y2, elem, n32, no_cyclic_elem, alpha, beta, gamma"""
  def __init__(self):
    self.path=''
    self.n2d=0
    self.n3d=0
    self.e2d=0
    self.nlev=0
    self.nodes=[]
    self.zlevs= []
    self.x2= []
    self.y2= []
    self.x3= []
    self.y3= []
    self.z3= []
    self.elem= []
    self.n32=[]
    self.no_cyclic_elem=[]
    self.topo=[]
    self.voltri=[]
    self.alpha=0
    self.beta=0
    self.gamma=0		
    def __str__(self):
      return "mesh path=%s" % self.path

def read_mesh(path, alpha, beta, gamma):
  import numpy as np
  from scalar_r2g import scalar_r2g

  mesh=fesom_mesh()
  mesh.path=path
  mesh.alpha=alpha
  mesh.beta=beta
  mesh.gamma=gamma	

  nod2dfile=mesh.path+'nod2d.out'
  elm2dfile=mesh.path+'elem2d.out'
  aux3dfile=mesh.path+'aux3d.out'
  nod3dfile=mesh.path+'nod3d.out'
  
  TmpArr = []
  x2,y2,nodes = [],[],[]
  reader = open(nod2dfile,'r')
  next (reader)
  for line in reader:
    TmpArr = line.split(' ')
    x2.append(float(TmpArr[1]))
    y2.append(float(TmpArr[2]))
    nodes.append(int(TmpArr[3])) 
  mesh.x2, mesh.y2, mesh.nodes, mesh.n2d = x2,y2,nodes,len(x2)	

  TmpArr = []
  x3,y3,z3,nodes = [],[],[],[]
  reader = open(nod3dfile,'r')
  next (reader)
  for line in reader:
    TmpArr = line.split(' ')
    if TmpArr[1].isdigit():
      x3.append(float(TmpArr[1]))
      y3.append(float(TmpArr[2]))
      z3.append(float(TmpArr[3]))
  mesh.x3, mesh.y3, mesh.z3, mesh.n3d, mesh.zlevs = x3,y3,z3,len(x3),np.unique(z3)

  TmpArr = []
  elem = []
  reader = open(elm2dfile,'r')
  for line in reader:
    TmpArr = line.split(' ')
    elem.append(int(TmpArr[0]))
  mesh.elem,mesh.e2d = elem, np.shape(mesh.elem)[0]

  with open(aux3dfile) as f:
    mesh.nlev=int(f.next())
    mesh.n32=np.array([f.next() for x in xrange(mesh.n2d*mesh.nlev)]).astype(int).reshape(mesh.n2d,mesh.nlev)	
    mesh.topo=np.zeros(shape=(mesh.n2d))
    for prof in mesh.n32:			
      ind_nan=[iz for iz in prof if iz>0]
      ind_nan=ind_nan[-1]
      mesh.topo[prof[0]-1]=z3[ind_nan-1]
    print type(x3)
    #we should rotate the mesh to the geographical coordinates
    (mesh.x3,mesh.y3)=scalar_r2g(alpha,beta,gamma,int(x3),int(y3))
    (mesh.x2,mesh.y2)=scalar_r2g(alpha,beta,gamma,x2,y2)
    d=mesh.x2[elem].max(axis=1)-mesh.x2[elem].min(axis=1)
    mesh.no_cyclic_elem = [i for (i, val) in enumerate(d) if val < 100]	

  return mesh
