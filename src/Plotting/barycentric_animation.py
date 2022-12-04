
def rot_func(coor):
    import numpy as np
    trans_mat = np.array([[1, np.cos(np.pi/3)],[0, np.sin(np.pi/3)]])
    ein_coor = np.array([coor[1],coor[0]])
    equil_coor = np.dot(trans_mat, ein_coor)
    return equil_coor

def barycentric(triangle, point):
    import numpy as np
    px = point[0];  py = point[1]
    x0 = triangle[0][0]; y0 = triangle[0][1]
    x1 = triangle[1][0]; y1 = triangle[1][1]
    x2 = triangle[2][0]; y2 = triangle[2][1]
    T = np.array([[x0-x2, x1-x2],[y0-y2, y1-y2]])
    Lambda = [0,0,0]
    Lambda[0] = ((y1-y2)*(px-x2) + (x2-x1)*(py-y2))/np.linalg.det(T)
    Lambda[1] = ((y2-y0)*(px-x2) + (x0-x2)*(py-y2))/np.linalg.det(T)
    Lambda[2] = 1 - Lambda[0] - Lambda[1]
    return Lambda

def transformedFromEisenstein(dual_faces, arcs, arcpos):
    data = []; data_trans = []
    import numpy as np
    for tri in dual_faces:
        tri_data = []; rot_tri = [];
        j = 0; tri_array = np.zeros([3,2]); b = 0
        while True:
            if tri[0] == arcs[j][0] and tri[1] == arcs[j][1]:
                tri_array[0,0] = rot_func(arcpos[j][0])[0];    tri_array[0,1] = rot_func(arcpos[j][0])[1]
                tri_array[1,0] = rot_func(arcpos[j][1])[0];    tri_array[1,1] = rot_func(arcpos[j][1])[1]
                
                k = 0
                while True:
                    if tri[1] == arcs[k][0] and tri[2] == arcs[k][1]:
                        tri_array[2,0] = rot_func(arcpos[k][1])[0];    tri_array[2,1] = rot_func(arcpos[k][1])[1]
                        b = 1
                        break
                    k+=1    
            j+=1
            if b == 1:
                break
        #data.append(tri_data);  #data_trans.append(rot_tri)
        data_trans.append(tri_array)
    return data_trans

#import BasisFunctions_Quadratic
#basis = BasisFunctions_Quadratic
#%%
def contineous_rep(data_trans, u, basis, grid_size_x, dual_faces):
    import numpy as np
    
    min_x = data_trans[0][0,0]; max_x = data_trans[0][0,0]
    min_y = data_trans[0][0,1]; max_y = data_trans[0][0,1]
    for i in range(len(data_trans)):
        for j in range(3):
            if data_trans[i][j,0] < min_x:
                min_x = data_trans[i][j,0]
            elif data_trans[i][j,0] > max_x:
                max_x = data_trans[i][j,0]
                
            if data_trans[i][j,1] < min_y:
                min_y = data_trans[i][j,1]
            elif data_trans[i][j,1] > max_y:
                max_y = data_trans[i][j,1]
    
    max_x *= 1.05; max_y *= 1.025
    init_grid_x = grid_size_x+10;
    init_grid_y = int(np.ceil((grid_size_x+10)*(max_y-min_y)/(max_x-min_x))) + 10
    #init_grid_y = np.sin(np.pi/3) * init_grid_x/1.6
    #init_grid_y = init_grid_x#
    #init_grid_y = grid_size_x * (np.abs(min_y)+np.abs(max_y))/(np.abs(min_x)+np.abs(max_x)) *np.sin(np.pi/3) 
    
    
    #init_grid_y = grid_size_x * (np.abs(min_y)+np.abs(max_y))/(np.abs(min_x)+np.abs(max_x)) + 10
    #init_grid_y = int(np.ceil(init_grid_y))  + 10
    
    #C20
    #init_grid_y = int(np.ceil(init_grid_x * 4/7 * np.sin(2/3 * np.pi)))
    
    #min_x = 1+np.cos(np.pi/3)*2; max_x = len(grid[0]) + np.cos(np.pi/3) * (len(grid) - 1)
    #min_y = np.sin(np.pi/3)
    pixel_spacing = (max_x - min_x)/init_grid_x
    
    #tri = data_trans[0]; tri_num = dual_faces[0]
    pixel_map = np.zeros([init_grid_y, init_grid_x])
    pixel_map -= 100
    #pixel_spacing = 6/init_grid_x;
    init = min_x
    for K in range(len(data_trans)):
        tri = data_trans[K];    tri_num = dual_faces[K]
        i_range = int(np.floor((np.max(tri[:,0]) - np.min(tri[:,0]))/pixel_spacing))
        j_range = int(np.floor((np.max(tri[:,1]) - np.min(tri[:,1]))/pixel_spacing))
        #for i in range(init_grid_x):
        start_x = int(np.floor((np.min(tri[:,0]) - min_x)/pixel_spacing))
        start_y = int(np.floor((np.min(tri[:,1]) - min_y)/pixel_spacing))
        #for i in range(start, start + i_range):
        for i in range(start_x-2, start_x + i_range+2):
            #for j in range(init_grid_y):
            for j in range(start_y - 2, start_y + j_range+2):
                mid_pixel = np.array([i*pixel_spacing + min_x, j*pixel_spacing + min_y])
                #mid_pixel = np.array([i*pixel_spacing + start, j*pixel_spacing + min_y])
                bary_var = barycentric(tri, mid_pixel)
                if bary_var[0] >= 0 and bary_var[0] <= 1:
                    if bary_var[1] >= 0 and bary_var[1] <= 1:
                        if bary_var[2] >= 0 and bary_var[2] <= 1:
                            #pixel_map[init_grid_y-j,i] = pixel_value(tri, bary_var, K, basis, u)
                            pixel_map[j,i] = pixel_value(tri, bary_var, K, basis, u, dual_faces)
        #print('Computing triangle ' + str (K))
    return pixel_map

#%%
def pixel_value(tri, bary_var, tri_num, basis, u, dual_faces):
    #import FEM_Assembly; 
    import numpy as np
    import FEM.FEM_Assembly as FEM_Assembly
    #dual_faces = get_input('dual_faces')
    Connectivity_matrix = basis.GlobalMapping(dual_faces)

    '''trans_mat = np.array([[tri[1,0] - tri[0,0], tri[2,0] - tri[0,0]], 
                          [tri[1,1] - tri[0,1], tri[2,1] - tri[0,1]]])
    inv_Jacobian = np.linalg.inv(trans_mat)
    quad_ref_sys = np.dot(inv_Jacobian,   coor - tri[0,:]) #Quadrature points in canonical system
    '''
    unit_tri = np.array([[0,0],[1,0],[0,1]])
    x = bary_var[0]*unit_tri[0,0] + bary_var[1]*unit_tri[1,0] + bary_var[2]*unit_tri[2,0]
    y = bary_var[0]*unit_tri[0,1] + bary_var[1]*unit_tri[1,1] + bary_var[2]*unit_tri[2,1]

    indices = FEM_Assembly.Local2Global(Connectivity_matrix, tri_num)
    val = 0
    for N in range(basis.LocalDOF()):
        #BasisFunctions_Quadratic.BasisFunctionValue(N, point)
        try:
            val += u[int(indices[N])] * basis.BasisFunctionValue(N, [x,y])
        except IndexError:
            print('What')
    return val

#%%
'''
dual_faces = get_input('dual_faces')
arcs = get_input('arcs')
arcpos = get_input('arcpos')
data_trans = transformedFromEisenstein(dual_faces, arcs, arcpos)
import matplotlib.pyplot as plt
plt.close('all')

fig, ax = plt.subplots()
for i in range(20):
    ax.scatter(data_trans[i][0,0], data_trans[i][0,1])
    ax.scatter(data_trans[i][1,0], data_trans[i][1,1])
    ax.scatter(data_trans[i][2,0], data_trans[i][2,1])
plt.show()

#%%
import numpy as np
u = np.zeros(42); import BasisFunctions_Quadratic
#u[4] = 1;   u[7] = 1;     u[5] = 1
u[9] = 1; u[10] = 1;
u[5] = 0.5
basis = BasisFunctions_Quadratic

from numpy.ma import masked_array

A = contineous_rep(data_trans, u, basis, 300)
A1 = masked_array(A, A == -100)
#A =np.copy(A1)
from matplotlib import cm
fig,ax = plt.subplots()

pb = ax.imshow(A1, interpolation='nearest',cmap=cm.inferno)
cba = plt.colorbar(pb,shrink=0.25)
'''
#%%
#%%

def get_input(types):
    if types == 'grid':
        grid = [[-1, 10, 10, 10, 10, 10],
                [2, 8, 9, 11, 1, 2],
                [3, 6, 7, 4, 0, 3],
                [5, 5, 5, 5, 5, -1]];        
        return grid
    elif types == 'arcs':
        arcs = [[2, 3], [3, 5], [6, 5], [7, 5], [4, 5], [0, 5], [5, 6], [5, 7], [5, 
  4], [5, 0], [5, 3], [8, 10], [9, 10], [11, 10], [1, 10], [2, 
  10], [3, 2], [10, 2], [10, 8], [10, 9], [10, 11], [10, 1], [2, 
  8], [8, 9], [9, 11], [11, 1], [1, 2], [3, 6], [6, 7], [7, 4], [4, 
  0], [0, 3], [8, 6], [9, 7], [11, 4], [1, 0], [3, 8], [6, 9], [7, 
  11], [4, 1], [0, 2], [8, 2], [9, 8], [11, 9], [1, 11], [2, 1], [6, 
  3], [7, 6], [4, 7], [0, 4], [3, 0], [6, 8], [7, 9], [4, 11], [0, 
  1], [8, 3], [9, 6], [11, 7], [1, 4], [2, 0]];
        return arcs

    elif types == 'arcpos':
        arcpos = [[[2, 1], [3, 1]], [[3, 1], [4, 1]], [[3, 2], [4, 2]], [[3, 3], [4, 
   3]], [[3, 4], [4, 4]], [[3, 5], [4, 5]], [[4, 1], [3, 2]], [[4, 
   2], [3, 3]], [[4, 3], [3, 4]], [[4, 4], [3, 5]], [[4, 5], [3, 
   6]], [[2, 2], [1, 2]], [[2, 3], [1, 3]], [[2, 4], [1, 4]], [[2, 
   5], [1, 5]], [[2, 6], [1, 6]], [[3, 6], [2, 6]], [[1, 2], [2, 
   1]], [[1, 3], [2, 2]], [[1, 4], [2, 3]], [[1, 5], [2, 4]], [[1, 
   6], [2, 5]], [[2, 1], [2, 2]], [[2, 2], [2, 3]], [[2, 3], [2, 
   4]], [[2, 4], [2, 5]], [[2, 5], [2, 6]], [[3, 1], [3, 2]], [[3, 
   2], [3, 3]], [[3, 3], [3, 4]], [[3, 4], [3, 5]], [[3, 5], [3, 
   6]], [[2, 2], [3, 2]], [[2, 3], [3, 3]], [[2, 4], [3, 4]], [[2, 
   5], [3, 5]], [[3, 1], [2, 2]], [[3, 2], [2, 3]], [[3, 3], [2, 
   4]], [[3, 4], [2, 5]], [[3, 5], [2, 6]], [[2, 2], [2, 1]], [[2, 
   3], [2, 2]], [[2, 4], [2, 3]], [[2, 5], [2, 4]], [[2, 6], [2, 
   5]], [[3, 2], [3, 1]], [[3, 3], [3, 2]], [[3, 4], [3, 3]], [[3, 
   5], [3, 4]], [[3, 6], [3, 5]], [[3, 2], [2, 2]], [[3, 3], [2, 
   3]], [[3, 4], [2, 4]], [[3, 5], [2, 5]], [[2, 2], [3, 1]], [[2, 
   3], [3, 2]], [[2, 4], [3, 3]], [[2, 5], [3, 4]], [[2, 6], [3, 5]]];
        return arcpos
    
    else:
        import numpy as np
        dual_faces = np.array([[6,7,9],[6,9,8],[5,7,6],[7,11,9],[4,7,5],[3,6,8],[2,8,10],[2,3,8],[9,11,10],[8,9,10],[4,11,7],[3,5,6],[1,10,11],[1,2,10],[1,11,4],[0,4,5],[0,1,4],[0,2,1],[0,3,2],[0,5,3]])
        return dual_faces
    return

