def fieldatdepth2d(depth,data3,mesh):
  import numpy as np
  d_ind  = np.nonzero((np.unique(mesh.z3) == depth))
  ep_ind  = np.squeeze(mesh.n32[:,d_ind])
  data2=np.zeros(shape=(mesh.n2d))
  ind_noempty=[i for i in range(mesh.n2d) if ep_ind[i]>=0]
  ind_empty=[i for i in range(mesh.n2d) if ep_ind[i]<0]
  data2[ind_noempty]=data3[ep_ind[ind_noempty]]
  data2[ind_empty]=np.nan
  return (data2)
