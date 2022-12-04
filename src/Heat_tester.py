import numpy as np
import os.path
folder = 'Coordinate_Quadrature_data'
coordinates_7 = np.load(os.path.join(folder, 'coordinates_54.npy'))
weights_7 = np.load(os.path.join(folder, 'weights_54.npy'))

Fullerene = 'C60nt'
index = 0
from FullereneData.c60nt import *
#%%
'''
dual_faces = np.array([[6,7,9],[6,9,8],[5,7,6],[7,11,9],[4,7,5],[3,6,8],[2,8,10],[2,3,8],[9,11,10],[8,9,10],[4,11,7],[3,5,6],[1,10,11],[1,2,10],[1,11,4],[0,4,5],[0,1,4],[0,2,1],[0,3,2],[0,5,3]])
grid = barycentric_animation.get_input('grid')
arcs = barycentric_animation.get_input('arcs')
arcpos = barycentric_animation.get_input('arcpos')
'''

dual_faces = (triangles[index]);    grid   = dual_neighbours[index]
arcs = unfolding_arcs[index];       arcpos = unfolding_arcpos[index]


#resolution = 100
#a = plot.plot_2D(rho, Basis, resolution, 'final with 10 orbitals')

#%%
#a = plot.plot_2D(orb[:,0]*np.conj(orb[:,0]), Basis, resolution)

#plot_orb(orb, 100)
#Implement stationary colormapping
def triangulation(pentagon_array, hexagon_array):
    import numpy as np
    if hexagon_array.shape[1] == 0:
        T = tri_C20(pentagon_array)
        return T
    N_T = 5*pentagon_array.shape[1] + 6*hexagon_array.shape[1]      #num of coloumns in triangulation
    T = np.zeros([3, N_T], dtype=int)                                          #Initialize array for the triangulation
    T_p = np.zeros([3,5*12])                                        #Triagnel array from pent
    T_h = np.zeros([3,6*hexagon_array.shape[1]])                    #Triagnel array from hex
        
    vertex_num = int(5/3*pentagons.shape[1] + 2*hexagons.shape[1])  #Current vertex num.
    vertex_num_max = vertex_num + 12 + hexagons.shape[1]            #Upper vertex num.
    iterater = 0
    while True:
        if iterater < 12:
            T_p[2,iterater*5:iterater*5+5] = vertex_num
            T_p[0,iterater*5:iterater*5+5] = pentagon_array[0][iterater]
            T_p[1,iterater*5:iterater*5+5] = np.roll(pentagon_array[0][iterater],1)

        T_h[2,iterater*6:iterater*6+6] = vertex_num + 12
        T_h[0,iterater*6:iterater*6+6] = hexagon_array[0][iterater]
        T_h[1,iterater*6:iterater*6+6] = np.roll(hexagon_array[0][iterater],1)
        
        vertex_num+=1;
        iterater+=1
        if vertex_num + 12 == vertex_num_max:
            break
    T[:,:5*12] = T_p;   T[:,5*12:] = T_h;
    return T

def tri_C20(pentagon_array):
    T_p = np.zeros([3,5*12])                                        #Triagnel array from pent
    vertex_num = int(5/3*pentagons.shape[1] + 2*hexagons.shape[1])  #Current vertex num.
    vertex_num_max = vertex_num + 12 + hexagons.shape[1]            #Upper vertex num.
    for iterater in range(12):
        T_p[2,iterater*5:iterater*5+5] = vertex_num
        T_p[0,iterater*5:iterater*5+5] = pentagon_array[0][iterater]
        T_p[1,iterater*5:iterater*5+5] = np.roll(pentagon_array[0][iterater],1)
        vertex_num+=1;
    return T_p


#%%

#dual_faces = np.array([[6,7,9],[6,9,8],[5,7,6],[7,11,9],[4,7,5],[3,6,8],[2,8,10],[2,3,8],[9,11,10],[8,9,10],[4,11,7],[3,5,6],[1,10,11],[1,2,10],[1,11,4],[0,4,5],[0,1,4],[0,2,1],[0,3,2],[0,5,3]])
'''
from numpy import array
pentagons    = array([[[0,1,2,3,4],[5,6,7,8,9],[10,11,12,13,14],[15,16,17,18,19],[20,21,22,23,24],[25,26,27,28,29],[30,31,32,33,34],[35,36,37,38,39],[40,41,42,43,44],[45,46,47,48,49],[50,51,52,53,54],[55,56,57,58,59]]]);

hexagons     = array([[[14,13,9,8,4,3],[19,18,10,14,3,2],[24,23,15,19,2,1],[29,28,20,24,1,0],[8,7,25,29,0,4],[12,34,33,5,9,13],[17,39,38,11,10,18],[22,44,43,16,15,23],[27,49,48,21,20,28],[6,54,53,26,25,7],[33,32,50,54,6,5],[11,38,37,30,34,12],[16,43,42,35,39,17],[21,48,47,40,44,22],[26,53,52,45,49,27],[31,56,55,51,50,32],[37,36,57,56,31,30],[42,41,58,57,36,35],[47,46,59,58,41,40],[52,51,55,59,46,45]]]);

pentagons    = array([[[0,1,2,3,4],[5,6,7,4,3],[8,9,5,3,2],[10,11,8,2,1],[12,13,10,1,0],[7,14,12,0,4],[15,16,14,7,6],[9,17,15,6,5],[11,18,17,9,8],[13,19,18,11,10],[14,16,19,13,12],[16,15,17,18,19]]]);

hexagons     = array([[]]);
T = (triangulation(pentagons, hexagons)).T
#%%
siz = Basis.GlobalDOF(T)
'''

#%%
def full_integration(v_h, quad_coor, quad_weights, triangulation, basis):
    from FEM.UnitCell_Computations import UnitCell
    full_integral = 0
    tri_glo = basis.GlobalMapping(triangulation)
    #tri_glo = triangulation
    for tri in tri_glo:
        vertex_values = v_h[tri]
        integral = 0
        for i in range(len(quad_weights)):
            val = 0
            point = UnitCell.CoordinateTransformation(quad_coor[i,:])
            for N in range(basis.LocalDOF()):
                val += v_h[tri[N]] * basis.BasisFunctionValue(N, point)
            integral += val * quad_weights[i]
        full_integral += integral
    return full_integral
#A = Basis.GlobalMapping(T)

#%%
import FEM.FEM_Assembly as FEM_Assembly
from scipy.sparse import linalg    #from scipy import linalg
#import  FEM.BasisFunctions.BasisFunctions_Quadratic as Basis
import  FEM.BasisFunctions.BasisFunctions_Linear as Basis
points = Basis.GlobalDOF(dual_faces)
#p = np.ones([42])*0.5; p[6] = 10;  p[7] = 10;  p[9] = 10; p[12] = 10;  p[13] = 10;  p[14] = 10
#p = np.zeros([42]); p[:] = 1
#p[10] = 1000
#p[:] = 1000

p = np.random.rand(points)
p[:] *= 100

print(full_integration(p, coordinates_7, weights_7, dual_faces, Basis))

#M, W, b = FEM_Assembly.Assemble(dual_faces, Basis, coordinates_7, weights_7, 'Heat', p, 'dual')
M, W, A_mat, b = FEM_Assembly.Assemble(dual_faces, Basis, coordinates_7, weights_7, 'Heat', p, 'dual')
M_array = M.toarray()
W_array = W.toarray()
A_mat_array = A_mat.toarray()
A_mat_inv = np.linalg.inv(A_mat_array)

triangulation = Basis.GlobalMapping(dual_faces)
u_old = np.ones([points]);
u_old = np.copy(p);#   b=0
print(u_old)
#%%
delta_time = 0.015; steps = 2000
for i in range(steps):        
    if i == 0:
        b = 0
    #print(np.min(u_old))
    #print(u_old[10])
    #lin_sys_vec = np.dot(M_array, u_old) - delta_time * np.dot(np.dot(np.linalg.inv(M_array), W_array), u_old) + b*delta_time
    lin_sys_vec = np.dot(M_array, u_old) - delta_time * np.dot(W_array, u_old) 
    #lin_sys_vec = M @ u_old - delta_time * W @ u_old + b
    u_new = np.linalg.lstsq(M_array, lin_sys_vec)[0]
    
    div = np.linalg.norm(u_old-u_new, 2)
    #print(full_integration(u_old, coordinates_7, weights_7, dual_faces, Basis))
    u_old = np.copy(u_new)  #Setting the heat source term to zero after the first step. 
print(u_old)


#%%
dual_faces = np.array([[6,7,9],[6,9,8],[5,7,6],[7,11,9],[4,7,5],[3,6,8],[2,8,10],[2,3,8],[9,11,10],[8,9,10],[4,11,7],[3,5,6],[1,10,11],[1,2,10],[1,11,4],[0,4,5],[0,1,4],[0,2,1],[0,3,2],[0,5,3]])
import FEM.FEM_Assembly as FEM_Assembly
from scipy.sparse import linalg    #from scipy import linalg
#import  FEM.BasisFunctions.BasisFunctions_Quadratic as Basis
import  FEM.BasisFunctions.BasisFunctions_Linear as Basis
points = Basis.GlobalDOF(dual_faces)
p = np.ones(points)
p[7] = 100
M_lin, W_lin, A_mat, b = FEM_Assembly.Assemble(dual_faces, Basis, coordinates_7, weights_7, 'Heat', p, 'dual')
M_lin = M_lin.toarray(); W_lin = W_lin.toarray()
u_lin = np.copy(b)

u_old = u_lin
u_lin_list = []

delta_time = 0.015; steps = 600
for i in range(steps):        
    if i == 0:
        b = 0
        
    if i == 0 or i == 20:
        u_lin_list.append(u_old)
    elif i == 70 or i == 500:
        u_lin_list.append(u_old)
    #lin_sys_vec = np.dot(M_array, u_old) - delta_time * np.dot(np.dot(np.linalg.inv(M_array), W_array), u_old) + b*delta_time
    #lin_sys_vec = np.dot(M_array, u_old) - delta_time * np.dot(W_array, u_old) 
    lin_sys_vec = np.dot(M_lin, u_old) - delta_time * np.dot(W_lin, u_old) 
    #lin_sys_vec = M @ u_old - delta_time * W @ u_old + b
    u_new = np.linalg.lstsq(M_lin, lin_sys_vec)[0]
    
    div = np.linalg.norm(u_old-u_new, 2)
    #print(full_integration(u_old, coordinates_7, weights_7, dual_faces, Basis))
    u_old = np.copy(u_new)  #Setting the heat source term to zero after the first step. 

#%% Torus Heat 
import FEM.FEM_Assembly as FEM_Assembly
import os.path
folder = 'Coordinate_Quadrature_data'
coordinates_7 = np.load(os.path.join(folder, 'coordinates_54.npy'))
weights_7 = np.load(os.path.join(folder, 'weights_54.npy'))

dual_faces = np.array([[0, 1, 9],  [0, 1, 2],  [0, 2, 73],  [0, 7, 73],  [0, 7, 72],  [0, 9, 72],  [1, 2, 3],  [1, 9, 10],  [1, 5, 10],  [1, 4, 5],  [1, 3, 4],  [2, 3, 26],  [2, 25, 26],  [2, 25, 73],  [3, 4, 79],  [3, 47, 79],  [3, 26, 47],  [4, 5, 6],  [4, 6, 8],  [4, 8, 79],  [5, 10, 78],  [5, 11, 78],  [5, 6, 11],  [6, 7, 11],  [6, 7, 8],  [7, 8, 73],  [7, 11, 72],  [8, 67, 73],  [8, 67, 79],  [9, 21, 72],  [9, 21, 37],  [9, 10, 37],  [10, 37, 46],  [10, 46, 78],  [11, 65, 78],  [11, 65, 72],  [12, 13, 15],  [12, 13, 14],  [12, 14, 21],  [12, 21, 22],  [12, 17, 22],  [12, 16, 17],  [12, 15, 16],  [13, 14, 36],   [13, 15, 77],   [13, 35, 77],  [13, 35, 36],  [14, 18, 36],  [14, 18, 37], [14, 21, 37], [15, 16, 64], [15, 64, 71], [15, 71, 77], [16, 17, 19], [16, 19, 20], [16, 20, 64], [17, 22, 66], [17, 23, 66], [17, 19, 23], [18, 20, 36], [18, 19, 20], [18, 19, 23], [18, 23, 37], [20, 36, 51], [20, 51, 64], [21, 22, 72], [22, 65, 72], [22, 65, 66], [23, 46, 66], [23, 37, 46], [24, 25, 74], [24, 25, 26], [24, 26, 38], [24, 28, 38], [24, 27, 28], [24, 27, 74], [25, 52, 74], [25, 52, 73], [26, 42, 47], [26, 38, 42], [27, 28, 29], [27, 29, 75], [27, 57, 75], [27, 57, 74], [28, 29, 39], [28, 38, 48], [28, 43, 48], [28, 39, 43], [29, 32, 39], [29, 30, 32], [29, 30, 75], [30, 31, 32], [30, 31, 76], [30, 58, 76], [30, 58, 75], [31, 32, 40], [31, 34, 40], [31, 33, 34], [31, 33, 76], [32, 39, 49], [32, 44, 49], [32, 40, 44], [33, 34, 35], [33, 35, 77], [33, 61, 77], [33, 61, 76], [34, 35, 41], [34, 40, 50], [34, 45, 50], [34, 41, 45], [35, 36, 41], [36, 41, 51], [38, 42, 80], [38, 48, 80], [39, 43, 81], [39, 49, 81], [40, 44, 82], [40, 50, 82], [41, 45, 83], [41, 51, 83], [42, 47, 53], [42, 53, 54], [42, 54, 80], [43, 48, 56], [43, 55, 56], [43, 55, 81], [44, 49, 59], [44, 59, 60], [44, 60, 82], [45, 50, 62], [45, 62, 63], [45, 63, 83], [46, 65, 66], [46, 65, 78], [47, 67, 79], [47, 53, 67], [48, 68, 80], [48, 56, 68], [49, 69, 81], [49, 59, 69], [50, 70, 82], [50, 62, 70], [51, 71, 83], [51, 64, 71], [52, 54, 74], [52, 53, 54], [52, 53, 67], [52, 67, 73], [54, 68, 74], [54, 68, 80], [55, 56, 57], [55, 57, 75], [55, 69, 75], [55, 69, 81], [56, 57, 68], [57, 68, 74], [58, 60, 76], [58, 59, 60], [58, 59, 69], [58, 69, 75], [60, 70, 76], [60, 70, 82], [61, 63, 77],[61, 62, 63], [61, 62, 70], [61, 70, 76], [63, 71, 77], [63, 71, 83]])
import  FEM.BasisFunctions.BasisFunctions_Linear as Basis
points = Basis.GlobalDOF(dual_faces)
#p = np.ones([42])*0.5; p[6] = 10;  p[7] = 10;  p[9] = 10; p[12] = 10;  p[13] = 10;  p[14] = 10
#p = np.zeros([42]); p[:] = 1
#p[10] = 1000
#p[:] = 1000

A = np.unique(dual_faces[np.where(dual_faces==0)[0]])
p = np.ones(points);
p[A] = 50
p[0] = 100

print(full_integration(p, coordinates_7, weights_7, dual_faces, Basis))

M, W, A_mat, b = FEM_Assembly.Assemble(dual_faces, Basis, coordinates_7, weights_7, 'Heat', p, 'dual')
M_array = M.toarray(); W_array = W.toarray()
A_mat_array = A_mat.toarray(); A_mat_inv = np.linalg.inv(A_mat_array)
W_array = np.dot(A_mat_inv, W_array)

triangulation = Basis.GlobalMapping(dual_faces)
u_old = np.ones([points]);
u_old = np.copy(p);#   b=0
#u_old = np.dot(M_array, u_old)
#%%
print(u_old)
delta_time = 0.015; steps = 2000
for i in range(steps):        
    if i == 0:
        b = 0
    if i == 10:
        print(u_old)
    elif i == 50:
        print(u_old)
    elif i == 100:
        print(u_old)
    elif i == 1000:
        print(u_old)
    lin_sys_vec = np.dot(M_array, u_old) - delta_time * np.dot(W_array, u_old) 
    u_new = np.linalg.lstsq(M_array, lin_sys_vec)[0]
    
    div = np.linalg.norm(u_old-u_new, 2)
    u_old = np.copy(u_new)  #Setting the heat source term to zero after the first step.

