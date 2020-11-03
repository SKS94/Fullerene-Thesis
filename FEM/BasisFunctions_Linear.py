'''
Linear basis functions defined on the right angled triangle
{0,0}{0,1}{1,0}. The function valued or the gradient can be obtained.

Defining more complex functions as basis then needs the same functions as below.
Used the quadratic functions for an example instead
'''

def GlobalDOF(dual_faces):
    import numpy as np
    #dof = dual_faces.shape[0]/2 + 2
    dof = np.max(dual_faces) + 1
    return int(dof)
def LocalDOF():
    return int(3)

def BasisFunctionGradient(N, point):
    import numpy as np
    diff_basis_xi_eta = np.array([[-1, -1], [1, 0], [0,1]])

    return diff_basis_xi_eta[N,:]

def BasisFunctionValue(N, point):
    if N==0:
        return 1-point[0]-point[1]
    elif N==1:
        return point[0]
    elif N==2:
        return point[1]

def GlobalMapping(dual_faces): 
    return dual_faces