'''
Module of functions that are do operations locally on the unit cell
'''
def LocalOverlap(quad_coor, quad_weights, BasisType, determinant):
    import numpy as np
    local_dof = BasisType.LocalDOF()
    M_local   = np.zeros([local_dof,local_dof])
    for i in range(local_dof):
        for j in range(local_dof):
            for q in range(quad_coor.shape[0]):
                point = CoordinateTransformation(quad_coor[q,:])
                #val += BasisType.BasisFunctionValue(i, point) * BasisType.BasisFunctionValue(j, point) *\
                    #quad_weights[q] / (np.sqrt(27)/4) * 1/2 * determinant
                M_local[i,j] += BasisType.BasisFunctionValue(i, point) * BasisType.BasisFunctionValue(j, point) *\
                                quad_weights[q] * determinant
                 #First the normalized to 1, then 1/2(unit triangle) then det(J)
    return M_local

def LocalGradientOverlap(quad_coor, quad_weights, Basis, determinant, Element):
    import numpy as np

    num_shape_func = Basis.LocalDOF()
    num_quad = len(quad_weights)
    W_local = np.zeros([num_shape_func, num_shape_func])

    for i in range(num_shape_func):
        for j in range(num_shape_func):
            for q in range(num_quad):
                diff_basis_i = Gradient(i, quad_coor[q,:], Basis, Element)
                diff_basis_j = Gradient(j, quad_coor[q,:], Basis, Element) #Different functions evaluated in the same quadrature point

                W_local[i,j] += np.dot(diff_basis_i, diff_basis_j) * quad_weights[q] * determinant #normalize to right angled, then determinant
    return W_local

def LocalPot(quad_coor, quad_weights, Basis, Input, determinant):
    import numpy as np; #import basis_type
    local_dof = Basis.LocalDOF()
    V_local_quad = np.zeros([local_dof,local_dof*len(quad_weights)])#i.e \phi_i(x_q)*phi_i(x_q) * w_q
    for w in range(len(quad_weights)):
        point = CoordinateTransformation(quad_coor[w,:])
        for i in range(local_dof):
            for j in range(local_dof):
                #pot_local_quad[i,j+local_dof*w] = Basis.BasisFunctionValue(i, point) * Basis.BasisFunctionValue(j, point) *\
                     #quad_weights[w] / (np.sqrt(27)/4) * 1/2 * det
                V_local_quad[i,j+local_dof*w] = Basis.BasisFunctionValue(i, point) * Basis.BasisFunctionValue(j, point) *\
                                                quad_weights[w] * determinant
    return V_local_quad

def LoadVectorMatrixUnit(quad_coor, quad_weights, Basis, determinant):
    import numpy as np; #import basis_type
    local_dof = Basis.LocalDOF()
    BasisWeightMatrix = np.zeros([local_dof,len(quad_weights)])#i.e \phi_i(x_q) * w_q
    for w in range(len(quad_weights)):
        point = CoordinateTransformation(quad_coor[w,:])
        for i in range(local_dof):
            #BasisWeightMatrix[i,w] = Basis.BasisFunctionValue(i,point) *\
                 #quad_weights[w] * (np.sqrt(27)/4) * 1/2 * det
            BasisWeightMatrix[i,w] = Basis.BasisFunctionValue(i,point) *\
                                     quad_weights[w] * determinant
    return BasisWeightMatrix

def CoordinateTransformation(variable):
    import numpy as np      #variable equal to coordinate in physical space or gradient in reference
    
    #from FEM.shape_construction import FE_construction
    from FEM.TriangularElements.shape_construction import FE_construction
    #det, Element = FE_construction('triangle_equil')
    Element = FE_construction('triangle_equil')[1]
    Jacobian = np.array([[Element[1,0] - Element[0,0],Element[2,0] - Element[0,0]], 
                        [Element[1,1] - Element[0,1], Element[2,1] - Element[0,1]]])

    inv_Jacobian = np.linalg.inv(Jacobian)
    quad_reference = np.dot(inv_Jacobian,   variable - Element[0,:]) #Quadrature points in canonical system
    return quad_reference

def Gradient(n, coor, Basis, Element):
    import numpy as np    
    quad_reference = CoordinateTransformation(coor) #Coordinates in reference system
    gradient_reference = Basis.BasisFunctionGradient(n, quad_reference)
    #gradient_physical = CoordinateAndGradientTransformation(gradient_reference, 'gradient', Element)

    Jacobian = np.array([[Element[1,0] - Element[0,0], Element[2,0] - Element[0,0]], 
                        [Element[1,1] - Element[0,1], Element[2,1] - Element[0,1]]])

    trans_Jacobian = np.transpose(Jacobian)
    inv_trans_Jacobian = np.linalg.inv(trans_Jacobian)
    gradient_physical = np.dot(inv_trans_Jacobian, gradient_reference)    
    return gradient_physical



'''
    for w in range(len(quad_weights)):
        point = CoordinateAndGradientTransformation(quad_coor[w,:], 'coordinate')
        for i in range(local_dof):
            for j in range(local_dof):
                pot_local_quad[i,j+local_dof*w] = BasisType.BasisFunctionValue(i, point) * BasisType.BasisFunctionValue(j, point) * quad_weights[w]
    return pot_local_quad
'''
'''
    for i in range(local_dof):
        for j in range(local_dof):
            val = 0
            for k in range(quad_coor.shape[0]):
                point = CoordinateAndGradientTransformation(quad_coor[k,:], 'coordinate')
                val += BasisType.BasisFunctionValue(i, point) * BasisType.BasisFunctionValue(j, point) * quad_weights[k]
            pot_local[i,j] = val
    return pot_local
'''
'''
def CoordinateAndGradientTransformation(variable, eval_type, *argv):
    import numpy as np      #variable equal to coordinate in physical space or gradient in reference
    if eval_type == 'coordinate':
        from FEM.shape_construction import FE_construction
        #det, Element = FE_construction('triangle_equil')
        Element = FE_construction('triangle_equil')[1]
        Jacobian = np.array([[Element[1,0] - Element[0,0],Element[2,0] - Element[0,0]], 
                          [Element[1,1] - Element[0,1], Element[2,1] - Element[0,1]]])

        inv_Jacobian = np.linalg.inv(Jacobian)
        quad_reference = np.dot(inv_Jacobian,   variable - Element[0,:]) #Quadrature points in canonical system
        return quad_reference

    elif eval_type == 'gradient':
        Element = argv[0]
        Jacobian = np.array([[Element[1,0] - Element[0,0], Element[2,0] - Element[0,0]], 
                          [Element[1,1] - Element[0,1], Element[2,1] - Element[0,1]]])

        gradient_reference = np.copy(variable)
        trans_Jacobian = np.transpose(Jacobian)
        inv_trans_Jacobian = np.linalg.inv(trans_Jacobian)
        gradient_physical = np.dot(inv_trans_Jacobian, gradient_reference)    
        return gradient_physical
'''
'''def Gradient(n, coor, Basis, Element):
    import numpy as np    
    quad_reference = CoordinateAndGradientTransformation(coor, 'coordinate') #Coordinates in reference system
    gradient_reference = Basis.BasisFunctionGradient(n, quad_reference)
    gradient_physical = CoordinateAndGradientTransformation(gradient_reference, 'gradient', Element)
    return gradient_physical'''

























#Not used I think, old implementation of the heat equation.
'''
def LoadVectorUnit(quad_coordinates, quad_weights, BasisType, K, K_numbers):
    import numpy as np;
    local_dof = BasisType.LocalDOF()
    load_vec_local = np.zeros([local_dof])
    for i in range(local_dof):
        for j in range(quad_coordinates.shape[0]):
            trans_coor = CoordinateAndGradientTransformation(quad_coordinates[j], 'coordinate')
            #if basis_type == 'quadratic':
                #load_vec_local[i] +=  basis_functions.shape_quadratic_order(i, trans_coor) * quad_weights[j] * rhs_function(K, quad_coordinates[j], K_numbers, quad_weights[j])
            load_vec_local[i] +=  BasisType.BasisFunctionValue(i, trans_coor) * quad_weights[j] * rhs_function(K, quad_coordinates[j], K_numbers, quad_weights[j])
            #else:
                #load_vec_local[i] +=  basis_functions.shape_linear_order(i, trans_coor) * quad_weights[j] * rhs_function(K, quad_coordinates[j], K_numbers, quad_weights[j])
            
    return load_vec_local
'''
'''
def rhs_function(tri_num, coor, K_numbers, w):
    import numpy as np
    if tri_num == 0 and (coor == np.array([0,0])).all():
        print('Function value at tri num 0 and x=0, y=0 is 10')
        return 1/w
    return 0
'''

'''def LocalGradientOverlap(quad_coor, quad_weights, Basis, det, Element):
    #import basis_functions
    import numpy as np

    num_shape_func = Basis.LocalDOF()
    num_quad = len(quad_weights)
    #T = np.ones([num_quad, num_shape_func**2])
    
    T_temp = np.zeros([num_shape_func, num_shape_func])
    #T_temp2 = np.zeros([num_shape_func, num_shape_func])
    for i in range(num_shape_func):
        #T_temp = np.ones([num_shape_func, num_shape_func])
        for j in range(num_shape_func):
            for q in range(num_quad):
                diff_basis_i = Gradient(i, quad_coor[q,:], Basis, Element)
                diff_basis_j = Gradient(j, quad_coor[q,:], Basis, Element) #Different functions evaluated in the same quadrature point
                #T_temp[i,j] += np.dot(diff_basis_i, diff_basis_j) * quad_weights[q]
                #T_temp[i,j] += np.dot(diff_basis_i, diff_basis_j) * quad_weights[q]/(np.sqrt(27)/4) *1/2 * np.sqrt(27)/2
                
                #T_temp[i,j] += np.dot(diff_basis_i, diff_basis_j) * quad_weights[q]/ \
                    #(np.sqrt(27)/4) * 1/2 * det     #normalize to right angled, then determinant

                T_temp[i,j] += np.dot(diff_basis_i, diff_basis_j) * quad_weights[q] * det     #normalize to right angled, then determinant

                #quad_ref_sys = CoordinateAndGradientTransformation(quad_coor[q,:], 'coordinate')
                #diff_basis_i =  Basis.BasisFunctionGradient(i, quad_ref_sys)
                #diff_basis_j =  Basis.BasisFunctionGradient(j, quad_ref_sys)
#                T_temp2[i,j] += np.dot(diff_basis_i, diff_basis_j) * quad_weights[q]/ \
                    #(np.sqrt(27)/4) * 1/2 * det     #normalize to right angled, then determinant
    return T_temp#, T_temp2
'''