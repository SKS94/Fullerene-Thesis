'''
Quadratic basis functions defined on the right angled triangle
{0,0}{0,1}{1,0}. The function valued or the gradient can be obtained.

Defining more complex functions as basis then needs the same functions as below.
Used the quadratic functions for an example instead
'''
def GlobalDOF(tri): #Returns global dof
    dof = (tri.shape[0]/2 + 2) + tri.shape[0]*3/2
    return int(dof)
def LocalDOF(): #Returns local dof
    return int(6)

def BasisFunctionGradient(N, point):    #Returns gradient in tranformed triangle
    import numpy as np
    if N==0:
        diff = np.array([-3 + 4*point[0] + 4*point[1],
                         -3 + 4*point[1] + 4*point[0]])
    elif N==2:
        diff = np.array([4*point[0] - 1, 0])
        
    elif N==4:
        diff = np.array([0, 4*point[1] - 1])
    
    elif N==1:
        diff = np.array([4 - 8*point[0] - 4*point[1], -4*point[0]])
    
    elif N==3:
        diff = np.array([4*point[1],    4*point[0]])
        
    elif N==5:
        diff = np.array([-4*point[1], 4 - 8*point[1] - 4*point[0]])
    return diff

def BasisFunctionValue(N, point): #Returns function value at transformed point
    if   N==0:
        return (1 - point[0] - point[1])*(1 - 2*point[0] - 2*point[1])
    elif N==2:
        return point[0]*(2*point[0] - 1)
    elif N==4:
        return point[1]*(2*point[1] - 1)
    elif N==1:
        return 4*point[0] * (1 - point[0] - point[1])
    elif N==3:
        return 4*point[0]*point[1]
    elif N==5:
        return 4*point[1] * (1 - point[0] - point[1])

'''def GlobalMapping(dual_faces): 
    
    #Creating a Connectivity matrix which extends the dual_faces to include the connecting
    #mid-points needed to describe the 6 quadratic functions on a unit cell
    
    import numpy as np
    ind_counter = np.max(dual_faces) + 1;    local_dof = LocalDOF()
    Connectivity_matrix = np.zeros([dual_faces.shape[0], local_dof], dtype=int)
    #Initializing the connectivity matrix with the given dual_faces
    Connectivity_matrix[:,0] = dual_faces[:,0]; Connectivity_matrix[:,2] = dual_faces[:,1]
    Connectivity_matrix[:,4] = dual_faces[:,2]
    #Initializing the first mid-points on cell
    Connectivity_matrix[0,1] = ind_counter; Connectivity_matrix[0,3] = ind_counter+1
    Connectivity_matrix[0,5] = ind_counter+2 

    ind_counter += 3;   pair = [0,1,2,0];  check = 0
    for k in range(1,dual_faces.shape[0]):     #Looping through each triangle
        for p in range(3):                 #Looping through the connecting corners in each triangle
            for ind in range(0, k):        #Looping to identify if the mid-points already has been assigned
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
    return Connectivity_matrix'''


def GlobalMapping(Tri): 
    import numpy as np
    ind_counter = np.max(Tri) + 1;    local_dof = LocalDOF()
    Connectivity_matrix = np.zeros([Tri.shape[0], local_dof], dtype=int)
    #Initializing the connectivity matrix with the given dual_faces
    Connectivity_matrix[:,0] = Tri[:,0]; Connectivity_matrix[:,2] = Tri[:,1]
    Connectivity_matrix[:,4] = Tri[:,2]
    #Initializing the first mid-points on cell
    Connectivity_matrix[0,1] = ind_counter; Connectivity_matrix[0,3] = ind_counter+1
    Connectivity_matrix[0,5] = ind_counter+2 

    ind_counter += 3;   pair = [0,1,2,0];  check = 0
    for k in range(1,Tri.shape[0]):               #Looping through each triangle
        for p in range(3):                 #Looping through the connecting corners in each triangle
            for ind in range(0, k):        #Looping to identify if the mid-points already has been assigned
                if Tri[k, pair[p]] in Tri[ind,:] and Tri[k, pair[p+1]] in Tri[ind,:]:
                    val1 = Tri[k, pair[p]];  val2 = Tri[k, pair[p+1]]
                        
                    ind1 = np.where(Tri[ind,:] == val1)[0]
                    ind2 = np.where(Tri[ind,:] == val2)[0]
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
