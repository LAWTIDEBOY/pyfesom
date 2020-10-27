def derivative(mesh):
  import numpy as np

  print 'Starting derivative calculation'
  elem_back = mesh.elem
  r_earth = 6367.5e3
  domain_length=360.
  domain_len=domain_length*np.pi/180.
  bafux_2d=np.zeros(shape=(np.shape(mesh.elem)))
  bafuy_2d=np.zeros(shape=(np.shape(mesh.elem)))
  el2d = np.shape(mesh.elem)
  el2d = el2d[0]
  voltriangle=np.zeros(shape=(el2d,1))
  full_cluster_area=np.zeros(shape=(el2d,1))

  derivative_stdbafu_x_2D= np.zeros(shape=(3,2))
  derivative_stdbafu_x_2D[0,:]= -1.
  derivative_stdbafu_x_2D[1,0]= 1.
  derivative_stdbafu_x_2D[2,1]= 1.

  local_cart = np.zeros(shape=(3,2))
  for i in range(0,el2d):
    node = elem_back[i,:]
    local_cart[:,0]=mesh.x3[node]*np.pi/180.  
    local_cart[:,1]=mesh.y3[node]*np.pi/180.

    # Jacobian
    jacobian2D = np.zeros(shape=(3,2)) 
    jacobian2D = local_cart[1:3,:]-np.vstack((local_cart[0,:],local_cart[0,:]))

    # check cyclic boundary
    for j in range(0,2):
      if (jacobian2D[j,0] > domain_len/3.0):
        jacobian2D[j,0] = jacobian2D[j,0]-domain_len
        if (jacobian2D[j,0] < -domain_len/3.0):
          jacobian2D[j,0] = jacobian2D[j,0]+domain_len

    jacobian2D = jacobian2D * r_earth
    jacobian2D[:,0] = jacobian2D[:,0] * np.mean(np.cos(local_cart[:,1]))

    #inverse
    jacobian2D_inv=np.linalg.inv(jacobian2D)
    derivative_locbafu_x_2D=np.dot(derivative_stdbafu_x_2D,np.matrix.transpose(jacobian2D_inv))
    bafux_2d[i,:] = derivative_locbafu_x_2D[:,0]
    bafuy_2d[i,:] = derivative_locbafu_x_2D[:,1]
    determinant   = np.linalg.det(jacobian2D)
    voltriangle[i] = np.abs(determinant)/2.0
    full_cluster_area[node] = full_cluster_area[node] + voltriangle[i]

  print 'Derivative calculated'   
  return(bafux_2d,bafuy_2d)
