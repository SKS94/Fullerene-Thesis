'''
README - This module can be called to construct the mass matrix, stiffness matrix and load vector 
of a fullerene given the dual faces. Pass the BasisType function set up as in the files BasisFunctions_Linear.py
and BasisFunctions_Linear.py. For a given function it might change whether or not all matrices are needed
so for the Poisson equation I will comment out the mass matrix to minimize calculations.
'''

'''
Rewritten version of stiffness_construction module.
The goal has been to make it more intuitive and non-repetative.
'''
def Assemble(Triangulation, Basis, quad_coordinates, quad_weights, equation_type, Input, *argv):
    '''The Assemble function which initialize all relevant assemblies for a given equation type.'''
    import numpy as np  
    quad_weights = quad_weights/(np.sqrt(27)/4) * 1/2

    global_dof = Basis.GlobalDOF(Triangulation);
    ConnectMatrix = Basis.GlobalMapping(Triangulation)
    
    from scipy.sparse.linalg import inv
    from scipy.sparse import lil_matrix

    if equation_type == 'Kohn-Sham':
        #M_global, W_global, Pot_eff = KohnShamSolver(Input, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, argv)
        #W_global = np.dot(A_Mat.toarray(), W_global.toarray())
        #M_global = np.dot(np.linalg.inv(A_Mat.toarray()), M_global.toarray())

        #A_mat_inv = inv(A_Mat)

        #Curved
        M_global, W_global, A_Mat, Pot_eff = KohnShamSolver(Input, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, argv)
        A_mat_inv = np.linalg.inv(A_Mat.toarray())
        W_global = A_mat_inv @ W_global
        W_global = lil_matrix((W_global))
        

        #W_global = A_Mat.dot(W_global)
        #M_global = A_mat_inv.dot(M_global)
        #W_global = A_Mat.dot(W_global)

        #W_global = A_mat_inv.dot(W_global)
        

        #Pot_eff = A_mat_inv.dot(Pot_eff)
        #Pot_eff =  np.dot(np.linalg.inv(A_Mat.toarray()), Pot_eff.toarray())

        return (1/2 * W_global + Pot_eff), M_global #A and B matrices

        #return (1/2 * W_global + Pot_eff), A_mat_inv #A and B matrices
        #return W_global, M_global, Pot_eff
            
    elif equation_type == 'Poisson':
        #P, f = PoissonSolver(Input*4*np.pi, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, argv)
        
        #W_global_cotan, f = PoissonSolver(Input*4*np.pi, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, argv)
        from scipy.sparse import lil_matrix
        
        #f = np.dot(np.linalg.inv(A_mat1), f)

        #curved
        
        W_global_cotan, A_Mat, f = PoissonSolver(Input*4*np.pi, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, argv)
        A_mat1 = A_Mat.toarray()
        f = A_Mat.dot(f)

        #f = np.dot(A_mat1, f)

        return W_global_cotan, f #P is the matrix and f is the load vector of the system
        #return P, f #P is the matrix and f is the load vector of the system

    elif equation_type == 'Heat':
        M, W, A, b = HeatSolver(Input, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, argv)
        #M, W, b = HeatSolver(Input, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, argv)
        return M, W, A, b
        #return M, W, b

    else:
        print('Equation name not valid')
        return


def KohnShamSolver(Input_pot, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, *argv):
    import numpy as np
    from scipy.sparse import lil_matrix
    M_global = lil_matrix((global_dof, global_dof))       #Mass matrix
    W_global = lil_matrix((global_dof, global_dof))       #Stiffness matrix
    V_global = lil_matrix((global_dof, global_dof))       #Stiffness matrix

    M_local_equil, M_local_isos =  MassLocal(Basis, quad_coordinates, quad_weights)    
    W_local_equil, W_local_isos = StiffnessLocal(Basis, quad_coordinates, quad_weights) #Constant for C20
    V_local_equil, V_local_isos = PotLocal(Basis, quad_coordinates, quad_weights, Input_pot)

    M_local = np.copy(M_local_isos)
    W_local = np.copy(W_local_isos)
    V_local = np.copy(V_local_isos)
    
    local_dof = Basis.LocalDOF()

    if 'dual' in argv[0]:
        V_local = np.copy(V_local_equil)
        W_local = np.copy(W_local_equil)
        M_local = np.copy(M_local_equil)
    for K in range(ConnectMatrix.shape[0]):
        if K == 12*5:
            V_local = np.copy(V_local_equil)
            W_local = np.copy(W_local_equil)
            M_local = np.copy(M_local_equil)

        glo_indices = Local2Global(ConnectMatrix, K)
        for ii in range(local_dof):
            for jj in range(local_dof):
                W_global._set_intXint(glo_indices[ii], glo_indices[jj], W_local[ii,jj] + W_global[glo_indices[ii], glo_indices[jj]])
                M_global._set_intXint(glo_indices[ii], glo_indices[jj], M_local[ii,jj] + M_global[glo_indices[ii], glo_indices[jj]])
                val = 0
                for w_q in range(len(quad_weights)): 
                    poly_expansion = contineous_rep(Input_pot, glo_indices, Basis, quad_coordinates[w_q,:]) #Polynomial expansion of the effective potential in a quadrature point
                    val += V_local[ii,jj + w_q*local_dof] * poly_expansion
                V_global._set_intXint(glo_indices[ii], glo_indices[jj], val + V_global[glo_indices[ii], glo_indices[jj]])
    W_global_cotan, A_mat = Cotangent_Laplacian(Basis, global_dof, ConnectMatrix)
    #W_global_cotan = np.dot(A_mat, W_global_cotan)

    return M_global, W_global_cotan, A_mat, V_global
    #return M_global, W_global_cotan, V_global

def HeatSolver(Input, BasisType, quad_coordinates, quad_weights, ConnectMatrix, global_dof, *argv):
    import numpy as np
    from scipy.sparse import lil_matrix

    local_dof = BasisType.LocalDOF()
    M_global = lil_matrix((global_dof, global_dof))       #Mass matrix
    W_global = lil_matrix((global_dof, global_dof))       #Stiffness matrix
    f_global = np.zeros([global_dof])

    M_local_equil, M_local_isos =  MassLocal(BasisType, quad_coordinates, quad_weights)    
    W_local_equil, W_local_isos = StiffnessLocal(BasisType, quad_coordinates, quad_weights) #Constant for C20
    Basis_w_mat_equil, Basis_w_mat_isos  = LoadVectorMatrix(quad_coordinates, quad_weights, BasisType)

    Basis_w_mat = np.copy(Basis_w_mat_isos)
    W_local = np.copy(W_local_isos)
    M_local = np.copy(M_local_isos)
    if 'dual' in argv[0]:
        Basis_w_mat = np.copy(Basis_w_mat_equil)
        W_local = np.copy(W_local_equil)
        M_local = np.copy(M_local_equil)
    for K in range(ConnectMatrix.shape[0]):
        if K == 12*5:
            Basis_w_mat = np.copy(Basis_w_mat_equil)
            W_local = np.copy(W_local_equil)
            M_local = np.copy(M_local_equil)

        glo_indices = Local2Global(ConnectMatrix, K)
        for ii in range(local_dof):
            for jj in range(local_dof):
                W_global._set_intXint(glo_indices[ii], glo_indices[jj], W_local[ii,jj] + W_global[glo_indices[ii], glo_indices[jj]])
                M_global._set_intXint(glo_indices[ii], glo_indices[jj], M_local[ii,jj] + M_global[glo_indices[ii], glo_indices[jj]])
            val = 0
            for w_q in range(len(quad_weights)):
                poly_expansion = contineous_rep(Input, glo_indices, BasisType, quad_coordinates[w_q,:]) #Polynomial expansion of the effective potential in a quadrature point
                val += Basis_w_mat[ii,w_q] * poly_expansion
            f_global[glo_indices[ii]] += val                
    W_global_cotan, A_mat = Cotangent_Laplacian(BasisType, global_dof, ConnectMatrix)
    return M_global, W_global_cotan, A_mat, f_global
    #return M_global, W_global, f_global

def PoissonSolver(Input_rho, BasisType, quad_coordinates, quad_weights, ConnectMatrix, global_dof, *argv):
    import numpy as np
    from scipy.sparse import lil_matrix
    local_dof = BasisType.LocalDOF()
    W_global = lil_matrix((global_dof, global_dof))       #Stiffness matrix
    W_local_equil, W_local_isos = StiffnessLocal(BasisType, quad_coordinates, quad_weights) #Constant for C20
    Basis_w_mat_equil, Basis_w_mat_isos = LoadVectorMatrix(quad_coordinates, quad_weights, BasisType)

    load_vec = np.zeros([global_dof])

    Basis_w_mat = np.copy(Basis_w_mat_isos)
    W_local = np.copy(W_local_isos)
    if 'dual' in argv[0]:
        Basis_w_mat = np.copy(Basis_w_mat_equil)
        W_local = np.copy(W_local_equil)

    for K in range(ConnectMatrix.shape[0]):
        if K == 12*5:
            Basis_w_mat = np.copy(Basis_w_mat_equil)
            W_local = np.copy(W_local_equil)

        glo_indices = Local2Global(ConnectMatrix, K)
        for ii in range(local_dof):
            for jj in range(local_dof):   
                W_global._set_intXint(glo_indices[ii], glo_indices[jj], W_local[ii,jj] + W_global[glo_indices[ii], glo_indices[jj]])
            val = 0
            for w_q in range(len(quad_weights)):
                poly_expansion = contineous_rep(Input_rho, glo_indices, BasisType, quad_coordinates[w_q,:]) #Polynomial expansion of the effective potential in a quadrature point
                val += Basis_w_mat[ii,w_q] * poly_expansion
            load_vec[glo_indices[ii]] += val                

    #return A_global, A_local, b_global, b_local
    W_global_cotan, Mat = Cotangent_Laplacian(BasisType, global_dof, ConnectMatrix)
    #return W_global, load_vec
    return W_global_cotan, Mat, load_vec

def contineous_rep(u, indices, Basis, point):
    local_dof = Basis.LocalDOF();   val_pol =  0
    import FEM.UnitCell_Computations.UnitCell as UnitCell
    point_transformed = UnitCell.CoordinateTransformation(point)
    for N in range(local_dof):
        val_pol += u[indices[N]] * Basis.BasisFunctionValue(N, point_transformed)
    return val_pol

def LoadVectorMatrix(quad_coordinates, quad_weights, Basis):
    #import UnitCell
    import FEM.UnitCell_Computations.UnitCell as UnitCell
    #import FEM.shape_construction as shape
    import FEM.TriangularElements.shape_construction as shape
    det, Element = shape.FE_construction('triangle_equil')
    LoadVecMat_equil = UnitCell.LoadVectorMatrixUnit(quad_coordinates, quad_weights, Basis, det)
    det, Element = shape.FE_construction('triangle_isosceles')
    LoadVecMat_isos = UnitCell.LoadVectorMatrixUnit(quad_coordinates, quad_weights, Basis, det)
    return LoadVecMat_equil, LoadVecMat_isos

def PotLocal(basis_type, quad_coordinates, quad_weights, Input):
    #import UnitCell
    import FEM.UnitCell_Computations.UnitCell as UnitCell
    #import FEM.shape_construction as shape
    import FEM.TriangularElements.shape_construction as shape
    det, Element = shape.FE_construction('triangle_equil')
    P_local_equil = UnitCell.LocalPot(quad_coordinates, quad_weights, basis_type, Input, det)
    det, Element = shape.FE_construction('triangle_isosceles')
    P_local_isos = UnitCell.LocalPot(quad_coordinates, quad_weights, basis_type, Input, det)
    return P_local_equil, P_local_isos

def Local2Global(Connectivity_matrix, tri_ind): #Local to global indices
    indices = Connectivity_matrix[tri_ind, :]
    return indices

def MassLocal(Basis, quad_coor, quad_weights):  #done
    import FEM.UnitCell_Computations.UnitCell as UnitCell
    #import FEM.shape_construction as shape
    import FEM.TriangularElements.shape_construction as shape
    det, Element = shape.FE_construction('triangle_equil')
    M_local_equil = UnitCell.LocalOverlap(quad_coor, quad_weights, Basis, det)
    det, Element = shape.FE_construction('triangle_isosceles')
    M_local_isos = UnitCell.LocalOverlap(quad_coor, quad_weights, Basis, det)
    return M_local_equil, M_local_isos 

def StiffnessLocal(Basis, quad_coor, quad_weights): #done
    import FEM.UnitCell_Computations.UnitCell as UnitCell
    #import FEM.shape_construction as shape
    import FEM.TriangularElements.shape_construction as shape
    #W_k = UnitCell.LocalGradientOverlap(quad_coor, quad_weights, BasisType)
    det, Element = shape.FE_construction('triangle_equil')
    W_local_equil = UnitCell.LocalGradientOverlap(quad_coor, quad_weights, Basis, det, Element)
    det, Element = shape.FE_construction('triangle_isosceles')
    W_local_isos = UnitCell.LocalGradientOverlap(quad_coor, quad_weights, Basis, det, Element)
    return W_local_equil, W_local_isos

def Cotangent_Laplacian(Basis, global_dof, dual_faces):
    '''Is the dual rep. set up where the first 12 indices determine pentagon faces?'''
    import numpy as np
    from scipy.sparse import lil_matrix
    Cotan_Laplacian = lil_matrix((global_dof, global_dof)) 
    A_mat = lil_matrix((global_dof, global_dof)) 

    Ab_pent = np.sqrt(27)/4 / 3 * 5 
    Ab_hex  = np.sqrt(27)/4 / 3 * 6 
    Ab_hep  = np.sqrt(27)/4 / 3 * 7 
    #Ab_pent = 5/3 * np.sqrt(3)/4 * (np.sqrt(3)/2)**2 / 4 #/4 new dual
    #Ab_hex  = 6/3 * np.sqrt(3)/4 * (np.sqrt(3)/2)**2 / 4 #/4 new dual
    w_ij = 1/np.sqrt(3)
    A = Ab_pent
    #A = 1
    for i in range(global_dof):
        #if i < 12:
        if len(np.where(dual_faces == i)[0]) == 5:
            A_mat._set_intXint(i, i,  Ab_pent)
        elif len(np.where(dual_faces == i)[0]) == 6:
            A_mat._set_intXint(i, i,  Ab_hex)
        elif len(np.where(dual_faces == i)[0]) == 7:
            A_mat._set_intXint(i, i,  Ab_hep)
        
        #if i == 12:
            #A = Ab_hex
            #A = 1
        triangle_ind = np.where(dual_faces == i)[0]
        ind = np.unique(dual_faces[triangle_ind])
        #Cotan_Laplacian[ind,i] = -1/np.sqrt(3) / A #(cot(2/6 * np.pi))
        for j in ind:
            Cotan_Laplacian._set_intXint(j, i,  - w_ij)
        '''
        if i < 12:
            #Cotan_Laplacian[i,i]  = 5/np.sqrt(3) / A
            Cotan_Laplacian._set_intXint(i, i,  5*w_ij)
        else:
            #Cotan_Laplacian[i,i]  = 6/np.sqrt(3) / A
            Cotan_Laplacian._set_intXint(i, i,  6*w_ij)
        '''
        if len(np.where(dual_faces == i)[0]) == 5:
            Cotan_Laplacian._set_intXint(i, i,  5*w_ij)
        elif len(np.where(dual_faces == i)[0]) == 6:
            Cotan_Laplacian._set_intXint(i, i,  6*w_ij)
        elif len(np.where(dual_faces == i)[0]) == 7:
            Cotan_Laplacian._set_intXint(i, i,  7*w_ij)
        

    
    return Cotan_Laplacian, A_mat




'''#Redundant, I think 
def LoadVectorLocal(quad_coordinates, quad_weights, K, K_numbers, Basis):
    #import UnitCell
    import FEM.UnitCell_Computations.UnitCell as UnitCell
    M_k = UnitCell.LoadVectorUnit(quad_coordinates, quad_weights, Basis, K, K_numbers)
    return M_k
'''
'''
def global_mapping(dual_faces, basis_type, global_dof):
    import numpy as np
    #if basis_type == 'linear':
        #return dual_faces
    elif basis_type == 'quadratic':
        ind_counter = np.max(dual_faces) + 1;    loc_dof = 6
        Connectivity_matrix = np.zeros([dual_faces.shape[0], loc_dof], dtype=int)
        Connectivity_matrix[:,0] = dual_faces[:,0]; Connectivity_matrix[:,2] = dual_faces[:,1]; Connectivity_matrix[:,4] = dual_faces[:,2]
        Connectivity_matrix[0,1] = ind_counter; Connectivity_matrix[0,3] = ind_counter+1; Connectivity_matrix[0,5] = ind_counter+2 #use ind_counter instead, #
        ind_counter += 3;   check = 0
        pair = [0,1,2,0]
        for k in range(1,dual_faces.shape[0]):  #Loops 
            for p in range(3):
                for ind in range(0, k):
                    if dual_faces[k, pair[p]] in dual_faces[ind,:] and dual_faces[k, pair[p+1]] in dual_faces[ind,:]:
                        val1 = dual_faces[k, pair[p]];  val2 = dual_faces[k, pair[p+1]]
                        
                        ind1 = np.where(dual_faces[ind,:] == val1)[0]
                        ind2 = np.where(dual_faces[ind,:] == val2)[0]
                        check = 1
                        if np.abs(ind1 - ind2) == 2:
                            Connectivity_matrix[k, pair[p]*2 +1] = Connectivity_matrix[ind, 5]
                        elif ind1 < ind2:
                            Connectivity_matrix[k, pair[p]*2 +1] = Connectivity_matrix[ind, ind1*2 + 1]
                        else:
                            Connectivity_matrix[k, pair[p]*2 +1] = Connectivity_matrix[ind, ind2*2 + 1]
                        break
                if check == 0:
                    Connectivity_matrix[k, p*2 + 1] = ind_counter
                    ind_counter += 1
                check = 0
    return Connectivity_matrix
'''
#%%

#import FEM_Assembly_sparse
'''
import BasisFunctions_Quadratic
import BasisFunctions_Linear
from shape_construction import FE_construction
import os.path

import numpy as np
#dual_faces = np.array([[6,7,9],[6,9,8],[5,7,6],[7,11,9],[4,7,5],[3,6,8],[2,8,10],[2,3,8],[9,11,10],[8,9,10],[4,11,7],[3,5,6],[1,10,11],[1,2,10],[1,11,4],[0,4,5],[0,1,4],[0,2,1],[0,3,2],[0,5,3]])
triangles_shape = [60,3];
#dual_faces = np.array([[[0,16,15],[0,15,14],[0,14,13],[0,13,12],[0,12,16],[1,17,22],[1,22,21],[1,21,16],[1,16,12],[1,12,17],[2,13,18],[2,18,23],[2,23,17],[2,17,12],[2,12,13],[3,14,19],[3,19,24],[3,24,18],[3,18,13],[3,13,14],[4,15,20],[4,20,25],[4,25,19],[4,19,14],[4,14,15],[5,16,21],[5,21,26],[5,26,20],[5,20,15],[5,15,16],[6,23,28],[6,28,27],[6,27,22],[6,22,17],[6,17,23],[7,24,29],[7,29,28],[7,28,23],[7,23,18],[7,18,24],[8,25,30],[8,30,29],[8,29,24],[8,24,19],[8,19,25],[9,26,31],[9,31,30],[9,30,25],[9,25,20],[9,20,26],[10,22,27],[10,27,31],[10,31,26],[10,26,21],[10,21,22],[11,31,27],[11,27,28],[11,28,29],[11,29,30],[11,30,31]]]).reshape(triangles_shape);
Element = FE_construction('triangle_equil')
folder = 'Coordinate_Quadrature_data'
coordinates_7 = np.load(os.path.join(folder, 'coordinates_7.npy'))
weights_7 = np.load(os.path.join(folder, 'weights_7.npy'))
'''
#%%
'''
W = Stiffness matrix,       M = Mass matrix,       b = Load vector
'''
'''W_sparse, b_sparse, M_sparse = Assemble(dual_faces, Element, BasisFunctions_Quadratic, coordinates_7, weights_7)
#W, b, M = FEM_Assembly.Assemble(dual_faces, Element, BasisFunctions_Linear, coordinates_7, weights_7)


#%%
u_old = np.zeros([b_sparse.shape[0]]);
u_old = np.copy(b_sparse.toarray()[:,0])
b = np.zeros([b_sparse.shape[0]])

delta_time = 0.05
steps = 200

for i in range(steps):        
    #lin_sys_vec = np.dot(M, u_old) - delta_time * np.dot(W, u_old) +  b
    
    lin_sys_vec = np.dot(M_sparse.toarray(), u_old) - delta_time * np.dot(W_sparse.toarray(), u_old) +  b
    #lin_sys_vec = M_sparse @ u_old - delta_time * (W_sparse @ u_old) +  b_sparse
    u_new = np.linalg.solve(M_sparse.toarray(), lin_sys_vec)
    
    
    div = np.abs(np.abs(u_old[0]) - np.abs(u_new[0]))
    if div < 0.000000001:
        break
    
    u_old = np.copy(u_new)  #Setting the heat source term to zero after the first step. 
    if i == 0:
        b = np.zeros([b_sparse.shape[0]])
        
    if i < 30:
        print(np.max(u_new))
    
'''

'''
def Assemble(Triangulation, Basis, quad_coordinates, quad_weights, equation_type, Input, *argv):
    #Constructs the matrices and ectors needed for the specfic equation.

    import numpy as np  
    quad_weights = quad_weights/(np.sqrt(27)/4) * 1/2

    global_dof = Basis.GlobalDOF(Triangulation);
    ConnectMatrix = Basis.GlobalMapping(Triangulation)
    
    #M_global = np.zeros([global_dof, global_dof])       #Mass matrix
    #W_global = np.zeros([global_dof, global_dof])       #Stiffness matrix
    #b_global = np.zeros([global_dof])                   #Load vector
    
    if equation_type == 'Kohn-Sham':
        M_global, W_global, Pot_eff = KohnShamSolver(Input, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, argv)
        return (1/2 * W_global + Pot_eff), M_global #A and B matrices
        #return W_global, M_global, Pot_eff

    elif equation_type == 'Poisson':
        P, f = PoissonSolver(Input*4*np.pi, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, argv)
        return P, f #P is the matrix and f is the load vector of the system

    elif equation_type == 'Heat':
        M, W, b = HeatSolver(Input, Basis, quad_coordinates, quad_weights, ConnectMatrix, global_dof, argv)
        return M, W, b

    else:
        print('Equation name not valid')
        return
    
    
    #nodes = ConnectMatrix[0,:] # Initialize local matrix with some nodes input
    W_local = StiffnessLocal(element, BasisType, quad_coordinates, quad_weights) #Constant for C20
    M_local =  MassLocal(element, BasisType, quad_coordinates, quad_weights)
    
    for K in range(tri_num):
        b_local = LoadVectorLocal(quad_coordinates, quad_weights, K, dual_faces[K,:], BasisType)
        
        nodes = ConnectMatrix[K,:]
        glo_indices = Local2Global(ConnectMatrix, K)
#        ind_count = 0;  test_ind = 6
        for ii in range(local_dof):
            #b_global[glo_indices[ii]] += b_local[ii]
            b_global._set_intXint(glo_indices[ii], 0, b_local[ii] + b_global[glo_indices[ii], 0])
            for jj in range(local_dof):
                #W_global[glo_indices[ii], glo_indices[jj]] += W_local[ii,jj]
                #M_global[glo_indices[ii], glo_indices[jj]] += M_local[ii,jj]
                
                W_global._set_intXint(glo_indices[ii], glo_indices[jj], W_local[ii,jj] + W_global[glo_indices[ii], glo_indices[jj]])
                M_global._set_intXint(glo_indices[ii], glo_indices[jj], M_local[ii,jj] + M_global[glo_indices[ii], glo_indices[jj]])
                

    #return A_global, A_local, b_global, b_local
    return W_global, b_global, M_global
    '''
    

