#First SCF-loop  for DFT

'''High negative repulsion to the effective, same sign as v_xc,
    Then when computed v_eff changes sign'''

#import BasisFunctions_Quadratic;    import BasisFunctions_Linear
#from shape_construction import FE_construction

def DFT_SCF(triangulation, tolerance, rho_initial, occupied_orb, Basis):
    fac = 1
    import Plotting.plot as plot;
    import numpy as np;    from scipy.sparse import linalg    #from scipy import linalg
    import FEM.FEM_Assembly as FEM_Assembly
    import os.path
    folder = 'Coordinate_Quadrature_data'

    coordinates_54 = np.load(os.path.join(folder, 'coordinates_54.npy'))#coordinates_7 = np.load(os.path.join(folder, 'coordinates_7.npy'))
    weights_54 = np.load(os.path.join(folder, 'weights_54.npy'))        #weights_7 = np.load(os.path.join(folder, 'weights_7.npy'))
    
    #weights_54 /= 4
    
    #triangulation = Basis.GlobalMapping(triangulation)
    counter = 0
    rho_const = full_integration(rho_initial, coordinates_54, weights_54, triangulation, Basis)
    print('Start rho constant is ' + str(rho_const))
    rho_old = rho_initial/rho_const * occupied_orb
    #rho_old = rho_old/rho_const * occupied_orb
    
    #rho_old /= np.sum(rho_old) 
    #rho_old *= 3
    '''All for limiting the hartree potential to integrate to same value'''
    '''
    c = 0.01
    v_h_sum, Mat = HartreeNormalization(c, occupied_orb, triangulation)
    v_h_sum = 1 * 50 / 2 
    W_modified = W_Poisson(coordinates_54, weights_54, Basis, rho_old, triangulation) #rho just used as input to FEM
    '''
    #W_modified = np.copy(Mat)
    #W_modified = W_Poisson(coordinates_54, weights_54, Basis, rho_old, triangulation, M)
    global_dof = Basis.GlobalDOF(triangulation);

    #full_integration(v_h, coordinates_54, weights_54, triangulation, Basis)
    #v_h_sum = 0.4618802153517056 * A_unit * triangulation.shape[0]
    
    all_nodes = Basis.GlobalMapping(triangulation)
    
    rho_full = []
    rho_full.append(rho_old)
    convergence = []
    energies_list = []
    while True:
        
        #a = plot.plot_2D(rho_old, Basis, resolution, 'Iteration num ' + str(counter), dual_faces, arcs, arcpos)
#        M, W, b = FEM_Assembly.Assemble(triangulation, Basis, coordinates_54, weights_54, 'Heat', rho_old, 'dual')
        #Solving the Poisson equation
        P, f = FEM_Assembly.Assemble(triangulation, Basis, coordinates_54, weights_54, 'Poisson', -rho_old, 'dual')
        #P1, f1 = FEM_Assembly.Assemble(triangulation, BasisFunctions_Quadratic, coordinates_7, weights_7, 'Poisson', rho_old)      
        #v_h = linalg.spsolve(P,f);        #v_h_alt = np.linalg.solve(P.toarray(), f);        #v_h = np.linalg.lstsq(P.toarray(), f)[0]
        
        #Mat = np.matrix([[P, np.identity(P.shape[0])], [M, 0]])        
        
        P, f = FEM_Assembly.Assemble(triangulation, Basis, coordinates_54, weights_54, 'Poisson', -rho_old, 'dual')
        least_square_solution = linalg.lsqr(P,f)[0]
        full = full_integration(least_square_solution, coordinates_54, weights_54, triangulation, Basis)
        f_0 = full/(np.sqrt(27)/4 * triangulation.shape[0])
        #f_0 = full/(np.sqrt(27)/4 * 20)
        #print(full_integration(least_square_solution - f_0, coordinates_54, weights_54, triangulation, Basis))
        v_h = least_square_solution
        v_h -= np.abs(np.max(v_h))
        #v_h = -least_square_solution
        
        #v_eff = v_effective(v_h, rho_old)
        
        v_eff = -v_effective(v_h, rho_old) #v_eff = v_effective(v_h, rho_old) INHERENTLY WRONG! DRAWN TOWARDS DENSITY
        #v_eff = v_effective(-v_h, rho_old) #USE   #Positive v_eff pushes away
        if counter == 20:
            print('stop')
        
        print('rho_new_energy = ' + str(full_integration(v_eff*rho_old, coordinates_54, weights_54, triangulation, Basis) +\
                         3/10*(3*np.pi)**(2/3)*full_integration(rho_old**(5/3), coordinates_54, weights_54, triangulation, Basis)))        
        E = full_integration(v_eff*rho_old, coordinates_54, weights_54, triangulation, Basis) +\
            3/10*(3*np.pi)**(2/3)*full_integration(rho_old**(5/3), coordinates_54, weights_54, triangulation, Basis)
        '''
        if counter != 0:
            E_T_ks = 0
            for q in range(occupied_orb):
                func = np.conj(orbitals_eigen_small[:,q]) * np.dot(P.toarray(), orbitals_eigen_small[:,q])
                E_T_ks += -1/2 * full_integration(func , coordinates_54, weights_54, triangulation, Basis)
                #print(-1/2 * full_integration(func , coordinates_54, weights_54, triangulation, Basis))
            #print('Energy in iteration ' + str(counter) + ' is ' + str(E_v_ks + E_T_ks))
            P
            #print('Energy in iteration ' + str(counter) + ' is ' + str(E_v_ks + E_T_ks))
            #print(E_T_rho)
            print('Energy in iteration alt' + str(counter) + ' is ' + str(E_v_ks + E_T_ks))
        '''

        #E_v_ks = full_integration(-v_eff*rho_old, coordinates_54, weights_54, triangulation, Basis)
        #E_T_rho = 3/10 *(6*np.pi**2)**(2/3) * full_integration(rho_old**(5/3), coordinates_54, weights_54, triangulation, Basis)
        #print('Energy in iteration ' + str(counter) + ' is ' + str(E_v_ks + E_T_rho))
        
        #Solving the Kohn-Sham equations
        #v_eff[:] = 0
        A, B = FEM_Assembly.Assemble(triangulation, Basis, coordinates_54, weights_54, 'Kohn-Sham', v_eff, 'dual')
        
#        B_inv = linalg.inv(B);#       B_inv_A = B_inv @ A
#        energy_eigen1, orbitals_eigen1 = linalg.eigs(B_inv_A, k=len(rho_initial)-2, which='SM') #Smallest magnitude is chose in solver
        
        energy_eigen, orbitals_eigen = linalg.eigs(A, k=len(rho_initial)-2, M=B, which='SM') #Smallest magnitude is chose in solver
        #energy_eigen, orbitals_eigen = linalg.eigs(A, k=200, M=B, which='SM') #Smallest magnitude is chose in solver
        
        index = np.argsort(energy_eigen)
        #Picking the eigenvectors/orbitals with the smallest energy/eigenvalue
        orbitals_eigen_small = orbitals_eigen[:, index[0:occupied_orb]]
        
        #Calculating the density
        orbitals_squared = (orbitals_eigen_small * np.conj(orbitals_eigen_small)).real
        #print(np.sum(energy_eigen[index[:occupied_orb]]))
        
        for k in range(occupied_orb):
            rho_norm = full_integration(orbitals_squared[:,k], coordinates_54, weights_54, triangulation, Basis)
            orbitals_squared[:,k] /=  rho_norm * fac
            #orbitals_squared[:,k] /=  10
            
            #print(full_integration(orbitals_squared[:,k], coordinates_54, weights_54, triangulation, Basis))
        rho_new = np.sum(orbitals_squared, axis=1)
        if counter > -1:
            alpha = 0.1
            if counter == 0:
                print('\n Density mixing initiallized \n')
            #if counter > 15: 
                #alpha = 0.15
            
            rho_new = alpha*rho_new + (1-alpha)*rho_old
        rho_norm = full_integration(rho_new, coordinates_54, weights_54, triangulation, Basis)
        rho_new = rho_new / rho_norm * occupied_orb
            
        '''
        if counter > 5:
            print('The norm of the difference is ' + str(np.linalg.norm(rho_new - rho_old, 2)))
            alpha = 0.5
            rho_new_alt = alpha*rho_new + (1-alpha)*rho_old
            rho_norm = full_integration(rho_new_alt, coordinates_54, weights_54, triangulation, Basis)
            rho_new_alt = rho_new_alt/rho_norm * occupied_orb
            rho_new = np.copy(rho_new_alt)
        '''
            
        #E_v_ks = full_integration(-v_eff*rho_new, coordinates_54, weights_54, triangulation, Basis)
        #E_T_rho = 3/10 *(6*np.pi**2)**(2/3) * full_integration(rho_new**(5/3), coordinates_54, weights_54, triangulation, Basis)
        #print('Energy in iteration ' + str(counter) + ' is ' + str(E_v_ks + E_T_rho))
        
        if counter == 45:
            v_x, v_c = pot_ExchangeCorrelation(rho_old)
            return rho_old, v_h, v_x, v_c, v_eff, rho_new
            
            #return rho_old, rho_new, orbitals_squared
            ind = [1,2,9,10,13,14]
            #orby = (orbitals_eigen[:,index[ind]] * np.conj(orbitals_eigen[:,index[ind]])).real
            energies =  energy_eigen[index[ind]].real
            return orbitals_eigen, energy_eigen, v_eff, index
            
            #return rho_old, v_h, v_x, v_c, v_eff, rho_new
        #rho_new /= np.sum(rho_new) 
        #rho_new *= 3
        
        '''rho_norm = full_integration(rho_new, coordinates_54, weights_54, triangulation, Basis)
        rho_new = rho_new/rho_norm * occupied_orb'''
        
        #rho_new =  rho_new/np.sum(rho_new)# /  occupied_orb #Normalize vector
        
        #print(np.sum(rho_new))
        
        
        rho_full.append(rho_new)
        convergence.append(np.linalg.norm(rho_new - rho_old, 2))
        energies_list.append(E)
        #energies_list.append(np.sum(energy_eigen[index[:occupied_orb]]))
        
        
        if counter == 10:
        #if counter == 5:
            #print('Maximum counter reached with norm ' + str(conver))
            #return rho_new, orbitals_squared
            energies =  energy_eigen[index[:occupied_orb]].real
            return rho_full, convergence, energies_list
            #return rho_new, convergence, energies
            #return rho_new, orbitals_squared, rho_old
            #return rho_new, convergence, energies_list
        
        
        if np.linalg.norm(rho_new - rho_old, 2) < tolerance:
            print(np.sum(energy_eigen[index[:occupied_orb]]))
            print('Convergence at ' + str(counter) + ' iterations')
            energies =  energy_eigen[index[:occupied_orb]].real
            #return rho_new, orbitals_squared, energies
            return rho_full, convergence, energies_list#, -v_h, v_eff
            #return rho_new, orbitals_squared, rho_old
            
            
            
            #return rho_new, convergence, energies_list
            #return rho_full, convergence
        else:
            conver = np.linalg.norm(rho_new - rho_old, 2)
            print('The norm of the difference is ' + str(np.linalg.norm(rho_new - rho_old, 2)))
            rho_old = np.copy(rho_new)
            #print('Tolerance not met with norm ' + str(np.sum(rho_new)))

        counter += 1
        '''if counter % 5 == 0:
            print('counter reached ' + str(counter))
            a = plot.plot_2D(rho_new, BasisFunctions_Quadratic, resolution)'''
    
#%% Testing the above function
def W_Poisson(quad_coor, quad_weights, Basis, rho, triangulation):
    import FEM.FEM_Assembly as FEM_Assembly
    H, M = FEM_Assembly.Assemble(triangulation, Basis, quad_coor, quad_weights, 'Kohn-Sham', rho, 'dual')
    P, f = FEM_Assembly.Assemble(triangulation, Basis, quad_coor, quad_weights, 'Poisson', rho, 'dual')
    global_dof = len(rho)
    N = P.shape[0]
    Mat = np.zeros([N*2, N*2])
    #Mat[:N, :N] = np.copy(P.toarray());   Mat[N:, :N] = np.copy(M.toarray()); 
    Mat[:N, :N] = np.copy(P.toarray());   Mat[N:, :N] = np.copy(M.toarray()); 
    Mat[:N, N:] = np.identity(N);   Mat[N:, N:] = np.zeros([N,N]); 
    return Mat

def HartreeNormalization(c, occupied_orb, triangulation):
    import os.path
    import numpy as np;    from scipy.sparse import linalg    #from scipy import linalg
    folder = 'Coordinate_Quadrature_data'
    coordinates_7 = np.load(os.path.join(folder, 'coordinates_7.npy'))
    weights_7 = np.load(os.path.join(folder, 'weights_7.npy'))
    
    coordinates_54 = np.load(os.path.join(folder, 'coordinates_54.npy'))
    weights_54 = np.load(os.path.join(folder, 'weights_54.npy'))
    import  FEM.BasisFunctions.BasisFunctions_Linear as Basis_init
    #import  FEM.BasisFunctions.BasisFunctions_Quadratic as Basis_init2
    import FEM.FEM_Assembly as FEM_Assembly
    A_unit = np.sqrt(3)/4 * 3
    rho_0_val = occupied_orb/(triangulation.shape[0]*A_unit);   
    rho_0 = np.ones([len(rho_initial)]); rho_0 *= rho_0_val
    H, M = FEM_Assembly.Assemble(triangulation, Basis_init, coordinates_54, weights_54, 'Kohn-Sham', rho_0, 'dual')
    P, f = FEM_Assembly.Assemble(triangulation, Basis_init, coordinates_54, weights_54, 'Poisson', -rho_0, 'dual')
    N = P.shape[0];     Mat = np.zeros([N*2, N*2])
    Mat[:N, :N] = np.copy(P.toarray());   Mat[N:, :N] = np.copy(M.toarray()); 
    Mat[:N, N:] = np.identity(N);   Mat[N:, N:] = np.zeros([N,N]); 
    vec = np.zeros(N * 2); vec[:N] = np.copy(f); vec[N:] = c
    least_square_solution = linalg.lsqr(P,f)
    v_h = np.linalg.lstsq(Mat,vec)[0][N:]
    v_h_sum = np.mean(v_h) * A_unit * triangulation.shape[0]
    return v_h_sum, Mat

def v_effective(v_h, rho):
    v_hartree = v_h
    v_x, v_c = pot_ExchangeCorrelation(rho)
    v_XC = v_x + v_c
    #v_XC = 0;     #v_XC = Exchange(rho) + Correlation(rho)      #v_XC = pot_ExchangeCorrelation(rho)
    v_ext = 0
    return (v_hartree + v_XC + v_ext)
    #return (v_hartree + v_XC + v_ext)

def pot_ExchangeCorrelation(rho):
    import numpy as np
    pot_XC = np.zeros(len(rho))
    v_x = Exchange(rho)
    v_c = Correlation(rho)
    pot_XC = v_x + v_c
    return v_x, v_c

def Exchange(rho):
    import numpy as np
    v_x = -(3/np.pi * rho)**(1/3)
    return v_x

def Correlation(rho):
    import numpy as np
    a = (np.log(2) - 1)/(2*np.pi**2)
    b = 20.4562557
    r_s = (3/(4*np.pi*rho))**(1/3)    
    c = (3/(4*np.pi))**(1/3)
    u = (1 + b*rho**(1/3)/c + b*rho**(2/3)/c**2) #By using the chain rule for the original expresseion of Chachiyo times density, recall that d/dx ln(u) = 1/u du/dx
    v_c = a*np.log(u) + a*rho*(2*b/(3*c**2*rho**(1/3)) + b/(3*c*rho**(2/3)))/u
#   return first_term + second_term
    return v_c

def full_integration(v_h, quad_coor, quad_weights, triangulation, basis):
    from FEM.UnitCell_Computations import UnitCell
    full_integral = 0
    tri_glo = basis.GlobalMapping(triangulation)
    tri_glo = tri_glo.astype(int)
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

def plot_orb(orb, res, Basis):
    global dual_faces, arcs, arcpos
    import Plotting.plot as plot
    for i in range(orb.shape[1]):
        #plot.plot_2D((orb[:,i]*np.conj(orb[:,i])).real, Basis, res, 'orbital num ' + str(i), dual_faces, arcs, arcpos)
        plot.plot_2D(orb[:,i], Basis, res, 'orbital num ' + str(i), dual_faces, arcs, arcpos)
    return

#%%
import numpy as np
'''
Fullerene = 'C120-D6'
index = 1
from FullereneData.c120_D6 import *

dual_faces = (triangles[index]);    grid   = dual_neighbours[index]
arcs = unfolding_arcs[index];       arcpos = unfolding_arcpos[index]
rho_initial = np.zeros([np.max(dual_faces)+1])
'''

'''
Fullerene = 'c60-nt';   index = 0
from FullereneData.c60nt import *

dual_faces = (triangles[index]);    grid   = dual_neighbours[index]
arcs = unfolding_arcs[index];       arcpos = unfolding_arcpos[index]
'''
import  FEM.BasisFunctions.BasisFunctions_Linear as Basis
#Torus
dual_faces = np.array([[0, 1, 9],  [0, 1, 2],  [0, 2, 73],  [0, 7, 73],  [0, 7, 72],  [0, 9, 72],  [1, 2, 3],  [1, 9, 10],  [1, 5, 10],  [1, 4, 5],  [1, 3, 4],  [2, 3, 26],  [2, 25, 26],  [2, 25, 73],  [3, 4, 79],  [3, 47, 79],  [3, 26, 47],  [4, 5, 6],  [4, 6, 8],  [4, 8, 79],  [5, 10, 78],  [5, 11, 78],  [5, 6, 11],  [6, 7, 11],  [6, 7, 8],  [7, 8, 73],  [7, 11, 72],  [8, 67, 73],  [8, 67, 79],  [9, 21, 72],  [9, 21, 37],  [9, 10, 37],  [10, 37, 46],  [10, 46, 78],  [11, 65, 78],  [11, 65, 72],  [12, 13, 15],  [12, 13, 14],  [12, 14, 21],  [12, 21, 22],  [12, 17, 22],  [12, 16, 17],  [12, 15, 16],  [13, 14, 36],   [13, 15, 77],   [13, 35, 77],  [13, 35, 36],  [14, 18, 36],  [14, 18, 37], [14, 21, 37], [15, 16, 64], [15, 64, 71], [15, 71, 77], [16, 17, 19], [16, 19, 20], [16, 20, 64], [17, 22, 66], [17, 23, 66], [17, 19, 23], [18, 20, 36], [18, 19, 20], [18, 19, 23], [18, 23, 37], [20, 36, 51], [20, 51, 64], [21, 22, 72], [22, 65, 72], [22, 65, 66], [23, 46, 66], [23, 37, 46], [24, 25, 74], [24, 25, 26], [24, 26, 38], [24, 28, 38], [24, 27, 28], [24, 27, 74], [25, 52, 74], [25, 52, 73], [26, 42, 47], [26, 38, 42], [27, 28, 29], [27, 29, 75], [27, 57, 75], [27, 57, 74], [28, 29, 39], [28, 38, 48], [28, 43, 48], [28, 39, 43], [29, 32, 39], [29, 30, 32], [29, 30, 75], [30, 31, 32], [30, 31, 76], [30, 58, 76], [30, 58, 75], [31, 32, 40], [31, 34, 40], [31, 33, 34], [31, 33, 76], [32, 39, 49], [32, 44, 49], [32, 40, 44], [33, 34, 35], [33, 35, 77], [33, 61, 77], [33, 61, 76], [34, 35, 41], [34, 40, 50], [34, 45, 50], [34, 41, 45], [35, 36, 41], [36, 41, 51], [38, 42, 80], [38, 48, 80], [39, 43, 81], [39, 49, 81], [40, 44, 82], [40, 50, 82], [41, 45, 83], [41, 51, 83], [42, 47, 53], [42, 53, 54], [42, 54, 80], [43, 48, 56], [43, 55, 56], [43, 55, 81], [44, 49, 59], [44, 59, 60], [44, 60, 82], [45, 50, 62], [45, 62, 63], [45, 63, 83], [46, 65, 66], [46, 65, 78], [47, 67, 79], [47, 53, 67], [48, 68, 80], [48, 56, 68], [49, 69, 81], [49, 59, 69], [50, 70, 82], [50, 62, 70], [51, 71, 83], [51, 64, 71], [52, 54, 74], [52, 53, 54], [52, 53, 67], [52, 67, 73], [54, 68, 74], [54, 68, 80], [55, 56, 57], [55, 57, 75], [55, 69, 75], [55, 69, 81], [56, 57, 68], [57, 68, 74], [58, 60, 76], [58, 59, 60], [58, 59, 69], [58, 69, 75], [60, 70, 76], [60, 70, 82], [61, 63, 77],[61, 62, 63], [61, 62, 70], [61, 70, 76], [63, 71, 77], [63, 71, 83]])

#rho_initial[23] = 1
#rho_initial += 0.00001

#rho_initial = np.zeros([242]); rho_initial += 0.01
#rho_initial[ConnectMatrix[0]] = 10

from untitled0 import Mesh_quadruple
import  FEM.BasisFunctions.BasisFunctions_Linear as Basis
'''
arcs = [[2, 3], [3, 5], [6, 5], [7, 5], [4, 5], [0, 5], [5, 6], [5, 7], [5, 
  4], [5, 0], [5, 3], [8, 10], [9, 10], [11, 10], [1, 10], [2, 
  10], [3, 2], [10, 2], [10, 8], [10, 9], [10, 11], [10, 1], [2, 
  8], [8, 9], [9, 11], [11, 1], [1, 2], [3, 6], [6, 7], [7, 4], [4, 
  0], [0, 3], [8, 6], [9, 7], [11, 4], [1, 0], [3, 8], [6, 9], [7, 
  11], [4, 1], [0, 2], [8, 2], [9, 8], [11, 9], [1, 11], [2, 1], [6, 
  3], [7, 6], [4, 7], [0, 4], [3, 0], [6, 8], [7, 9], [4, 11], [0, 
  1], [8, 3], [9, 6], [11, 7], [1, 4], [2, 0]];
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
dual_faces = np.array([[6,7,9],[6,9,8],[5,7,6],[7,11,9],[4,7,5],[3,6,8],[2,8,10],[2,3,8],[9,11,10],[8,9,10],[4,11,7],[3,5,6],[1,10,11],[1,2,10],[1,11,4],[0,4,5],[0,1,4],[0,2,1],[0,3,2],[0,5,3]])
'''
ConnectMatrix = Basis.GlobalMapping(dual_faces)
#dual_faces, arcs, arcpos = Mesh_quadruple(dual_faces, arcs, arcpos)
#dual_faces, arcs, arcpos = Mesh_quadruple(dual_faces, arcs, arcpos)
#dual_faces, arcs, arcpos = Mesh_quadruple(dual_faces, arcs, arcpos)
#dual_faces, arcs, arcpos = Mesh_quadruple(dual_faces, arcs, arcpos)

rho = np.random.rand(Basis.GlobalDOF(dual_faces))
#dual_faces_new, arcs_new, arcpos_new = Mesh_quadruple(dual_faces, arcs, arcpos)
import  FEM.BasisFunctions.BasisFunctions_Linear as Basis
#import  FEM.BasisFunctions.BasisFunctions_Quadratic as Basis
#rho_initial = np.zeros(Basis.GlobalDOF(dual_faces)) + 0.01
#rho_initial = np.ones(Basis.GlobalDOF(dual_faces))
rho_initial = np.random.rand(Basis.GlobalDOF(dual_faces))
#rho_initial[0] = 10;
#rho_initial[5] = 10;
#rho_initial[[26,27,29,30,33]] = 5
#%%
'''
for i in range(42, 162):
    ind =  np.where(T==i)[0]
    for j in range(2):
        tri = T[ind[j]]
        spec_ind = np.where(tri == i)[0]
        val = np.abs(rho_initial[np.roll(tri, 1)[spec_ind]] + rho_initial[np.roll(tri, -1)[spec_ind]])
        rho_initial[i] = val/2
#rho_initial[[17,18,22,32,38]] = 5
#rho_initial[[12,7,13, 42, 43, 44]] = 100
#rho_initial[[12,7,13]] = 100
#rho_initial[[6,7,9]] = 10
#rho_initial[[42,7,43]] = 100
#rho_initial[:12] *= 6/5
'''
#%%
T = Basis.GlobalMapping(dual_faces)
import Plotting.plot as plot
#rho, orb = DFT_SCF(T, 0.0001, rho_initial, 30)
orb_num = 20
#10 on C60 worked(one refinement)
resolution = 50

'''
full_rho = [];  full_conver = [];   full_energy = []
for i in range(1, 26):
#rho, orb = DFT_SCF(dual_faces, 0.00000001, rho_initial, orb_num, Basis)
    print(i)
    rho, conver_list, energy_list = DFT_SCF(dual_faces, 0.000001, rho_initial, i, Basis)
    full_rho.append(rho)
    full_conver.append(conver_list)
    full_energy.append(energy_list)
'''

'''remember number of solutions
https://www.diva-portal.org/smash/get/diva2:864857/FULLTEXT01.pdf'''

#plot.plot_2D(v_h, Basis, resolution, 'v_h ' + str(orb_num) +  ' Orbitals', dual_faces, arcs, arcpos)
#plot.plot_2D(v_eff, Basis, resolution, 'V_eff ' + str(orb_num) +  ' Orbitals', dual_faces, arcs, arcpos)
#%%
#rho, orb, rho_old, v_h, v_eff = DFT_SCF(dual_faces, 0.0001, rho_initial, orb_num, Basis)
#rho_list_120_1, convergence_120_1, energy_120_1 = DFT_SCF(dual_faces, 0.000001, rho_initial, 70, Basis)
rho_list_120_2, convergence_120_2, energy_120_2 = DFT_SCF(dual_faces, 0.00025, rho_initial, 1, Basis)
#%%
#plot.plot_2D(rho[-1], Basis, resolution, 'v_h ' + str(orb_num) +  ' Orbitals', dual_faces, arcs, arcpos)
#%%
#rho_list, convergence, energy = DFT_SCF(dual_faces, 0.01, rho_initial, 1, Basis)
#Torus
rho = []; con = []; en = []
#rho_torus, convergence_torus, energy_torus = DFT_SCF(dual_faces, 0.001, rho_initial, 1, Basis)
for i in range(40):
    rho_torus, convergence_torus, energy_torus = DFT_SCF(dual_faces, 0.000001, rho_initial, i+1, Basis)
    rho.append(rho_torus); con.append(convergence_torus); en.append(energy_torus)
#%%
rho_full_list = []; convergence_full_list = []
energy_full_list = []
for i in range(1,31):
    rho_list, convergence, energy = DFT_SCF(dual_faces, 0.0001, rho_initial, i, Basis)
    rho_full_list.append(rho_list)
    convergence_full_list.append(convergence)
    energy_full_list.append(energy)
    
    print(' \n ' + str(i) + ' reached \n')
my_list = [rho_full_list, convergence_full_list, energy_full_list]
with open('SCF_C20_data.txt', 'w') as f:
    f.write('First all rho, convergence then energy')
    for item in my_list:
        f.write('\n New data\n ')
        f.write("%s\n" % item)

#%%
plt.figure()
for j in range(15):
    x = len(convergence_full_list[j])
    x = np.linspace(1,x,x)
    plt.plot(x, convergence_full_list[j], label = j)
plt.yscale('log')
plt.legend(loc="upper right", ncol = 3, fontsize = 15)
plt.ylim([10**(-2),15])

plt.figure()
for j in range(15,30):
    x = len(convergence_full_list[j])
    x = np.linspace(1,x,x)
    plt.plot(x, convergence_full_list[j], label = j)
plt.yscale('log')
plt.legend(loc="upper right", ncol = 3, fontsize = 15)
plt.ylim([10**(-2),15])
    #%%
for j in range(30):
    if len(energy_full_list[j]) != 42:
        plt.figure()
        u = rho_full_list[j][-1]
        data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
        A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
        A1 = masked_array(A, A == -100)
        pb = plt.imshow(A1, interpolation='nearest', cmap=cm.inferno, origin='lower', vmin = 0)
        #ax_lst[i].set_axis_off()
        plt.colorbar(pb)
        plt.axis('off')
        plt.title('DFT Simulation using orbitlas num ' + str(j+1))
        
'''Show results for 4, 8, 25, 
look at 5 for a non-converged'''        
#%%
plot.plot_2D(rho_list[0], Basis, 300, 'Final Density Consisting of ' + str(orb_num) +  ' Orbitals', dual_faces, arcs, arcpos)
plot.plot_2D(rho_list[1], Basis, 300, 'Final Density Consisting of ' + str(orb_num) +  ' Orbitals', dual_faces, arcs, arcpos)
plot.plot_2D(rho_list[-1], Basis, 300, 'Final Density Consisting of ' + str(orb_num) +  ' Orbitals', dual_faces, arcs, arcpos)


#%%
plot.plot_2D(rho, Basis, 300, 'Final Density Consisting of ' + str(orb_num) +  ' Orbitals', dual_faces, arcs, arcpos)
plot.plot_2D(rho_old, Basis, 300, 'Final Density Consisting of ' + str(orb_num) +  ' Orbitals', dual_faces, arcs, arcpos)

plot_orb(orb[:,:orb_num], 300, Basis)

plot.plot_2D(rho, Basis, 300, 'Final Density Consisting of ' + str(orb_num) +  ' Orbitals', dual_faces, arcs, arcpos)
#rho, orb = DFT_SCF(dual_faces, 0.0000001, rho_initial, orb_num, Basis) #Quad C20

#rho_f, con = DFT_SCF(dual_faces, 0.001, rho_initial, orb_num, Basis) #Quad C20

#rho_old, v_h, v_x, v_c, v_eff, rho_new = DFT_SCF(dual_faces, 0.001, rho_initial, orb_num, Basis) #Quad C20
#plot_pic(rho_old, v_h, v_x, v_c, v_eff, rho_new, Basis, dual_faces, arcs, arcpos, 100)

#rho_old, rho_new, orb = DFT_SCF(dual_faces, 0.001, rho_initial, orb_num, Basis) #Quad C20
#orbs, energies, v, indices = DFT_SCF(dual_faces, 0.001, rho_initial, orb_num, Basis) #Quad C20
#%%
import Plotting.plot as plot
import time


#import  FEM.BasisFunctions.BasisFunctions_Quadratic as Basis
#import  FEM.BasisFunctions.BasisFunctions_Linear as Basis
#print(rho)
orb_num_list = [70, 100, 130]
for i in range(1):
    #resolution = 10*(i+1)
    #orb_num = 100
    orb_num = orb_num_list[i]
    rho, orb, energies = DFT_SCF(dual_faces, 0.00001, rho_initial, orb_num, Basis)
    
    resolution = 400
    t0 = time.time()
    a = plot.plot_2D(rho, Basis, resolution, 'Final Density Consisting of ' + str(orb_num) +  ' Orbitals', dual_faces, arcs, arcpos)
    t1 = time.time()

    total = t1-t0
    print(total)
#%%
#plot.plot_2D(orb[:,:5], Basis, resolution,'')
#orb = (orbs[:,indices]*np.conj(orbs[:,indices])).real
plot_orb(orb[:,:orb_num], 50, Basis)
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
'''
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
#T = (triangulation(pentagons, hexagons)).T
#siz = Basis.GlobalDOF(T)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
from math import sin, cos, pi
def regularPoly(n,a,b,r):
    points = [(a,b+r)]
    theta = pi/2
    dTheta = 2*pi/n
    for i in range(1,n):
        theta += dTheta
        points.append((a + r*cos(theta), b + r*sin(theta)))
    return points

P = regularPoly(5, 0, 0, 0.8506508083520399)

H = regularPoly(6, 0, 0, 1)

x = [P[0][0], P[1][0], P[2][0], P[3][0], P[4][0]]
y = [P[0][1], P[1][1], P[2][1], P[3][1], P[4][1]]
z = [0,0,0,0,0]
verts = [list(zip(x,y,z))]
ax.add_collection3d(Poly3DCollection(verts, alpha = 0.5))

Px = np.mean(x)
Py = np.mean(y)
#angle = np.pi*2/6
angle = np.pi*5/24
Rx = np.array([[1,0,0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])
angle_z = 2*np.pi/12
Rz = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], [np.sin(angle_z), np.cos(angle_z), 0], [0,0,1]])
#for i in range(5):
    #ax.scatter(P[i][0], P[i][1], 0, c='r')

cc = []
for j in range(6):
    coor = [H[j][0], H[j][1], 0]
    coor = np.dot(Rz, coor)
    c = np.dot(Rx, coor)
    if j == 0:
        val1 = c[2]
        val2 = c[1]
    print(c)
    c[2] -=   val1
    c[1] -= 0.6881909602355867 + val2
    #print(c)
    cc.append(c)
    #ax.scatter(c[0], c[1], c[2], c='b')
   
dualx = []; dualy = []; dualz = []

for k in range(2,7):
    angle_z = 2*np.pi/5 * k
    Rz = np.array([[np.cos(angle_z), -np.sin(angle_z), 0], [np.sin(angle_z), np.cos(angle_z), 0], [0,0,1]])
    x = []; y = []; z = [];
    for i in range(6):
        c1 = np.dot(Rz,cc[i])
        ax.scatter(c1[0], c1[1], c1[2], c='g', alpha = 0.0)
        x.append(c1[0]);    y.append(c1[1]);    z.append(c1[2]);
    verts = [list(zip(x,y,z))]
    ax.add_collection3d(Poly3DCollection(verts, alpha = 0.2, edgecolors='black', facecolor='orange'))
    ax.scatter(np.mean(x), np.mean(y), np.mean(z), c='orange', s= 40)
    dualx.append(np.mean(x));   dualy.append(np.mean(y));   dualz.append(np.mean(z))

ax.scatter(Px,Py,0,c='blue', s=40)
for i in range(5):
    ax.plot(dualx[i:i+2], dualy[i:i+2], dualz[i:i+2], c='black', linestyle='dashed')
    ax.plot([dualx[i], Px], [dualy[i], Py], [dualz[i], 0], c='black', linestyle='dashed')
ax.plot([dualx[-1], dualx[0]], [dualy[-1], dualy[0]], [dualz[-1], dualz[0]], c='black', linestyle='dashed')

plt.axis('off')
ax.set_xlim([-2.5,2.5])    
'''  
#%%
def plot_pic(u0, u1, u2, u3, u4, u5, Basis, dual_faces, arcs, arcpos, resolution):
    from Plotting import barycentric_animation
    from numpy.ma import masked_array
    from matplotlib import pyplot as plt
    from matplotlib import cm
    import matplotlib
    fig = plt.figure()
    #fig.suptitle('Heat dissipation in C20')
    plt.subplot(221)
    #norm1 = matplotlib.colors.Normalize(vmin=20, vmax=100)
    for i in range(6):
        ax = plt.subplot(3,2,i+1)
        if i == 0:
            u = np.copy(u0)
            ax.set_title(r'(a) $\rho_{old}$')
        elif i == 1:
            u = np.copy(u1)
            ax.set_title(r'(b) $v_H$')
        elif i == 2:
            u = np.copy(u2)
            ax.set_title(r'(c) $v_x$')
        elif i == 3:
            u = np.copy(u3)
            ax.set_title(r'(d) $v_c$')
        elif i == 4:
            u = np.copy(u4)
            ax.set_title(r'(e) $v_{KS}$')
        elif i == 5:
            u = np.copy(u5)
            ax.set_title(r'(f) $\rho_{new}$')
        data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
        A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
        A1 = masked_array(A, A == -100)
    #A =np.copy(A1)
    
    #fig,ax = plt.subplots()
        pb = ax.imshow(A1, interpolation='nearest', cmap=cm.inferno)
        #fig.colorbar(pb, cax=cbar, orientation='vertical')
        plt.colorbar(pb)
        plt.axis('off')
        #
    #fig.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #fig.colorbar(pb, cax=cbar_ax)
    plt.show()
    return
#%%
plot_pic(rho_old, v_h, v_x, v_c, v_eff, rho_new, Basis, dual_faces, arcs, arcpos, 100)


#%%
def plot_pic_alt(Basis, dual_faces, arcs, arcpos, resolution, Basis1):
    from Plotting import barycentric_animation;     from numpy.ma import masked_array
    from matplotlib import pyplot as plt;     from matplotlib import cm
    import matplotlib
    rho_initial = np.zeros([42])
    rho_initial[:] =  0; rho_initial += 0.5
    rho_initial[[5, 17, 7, 12, 6, 18]] = 10
    u0, u1, v_x, v_c, v_eff, rho_new = DFT_SCF(dual_faces, 0.001, rho_initial, orb_num, Basis) #Quad C20
    rho_initial = np.array([0.12344983, 0.59663442, 0.44264023, 0.92279873, 0.84962093,
       0.63840221, 0.69697079, 0.62450609, 0.58068021, 0.97640794,
       0.57323262, 0.24179967, 0.6479123 , 0.91338887, 0.13374419,
       0.71672089, 0.97985304, 0.02642397, 0.33265889, 0.07181509,
       0.35524093, 0.55692565, 0.20448581, 0.76462363, 0.23695562,
       0.8079494 , 0.97141425, 0.1513139 , 0.06386702, 0.17184362,
       0.47831249, 0.72574467, 0.2201492 , 0.74890521, 0.13521605,
       0.88560404, 0.72664812, 0.56765007, 0.44514012, 0.97878308,
       0.91846121, 0.95290903])
    u2, u3, v_x, v_c, v_eff, rho_new = DFT_SCF(dual_faces, 0.001, rho_initial, orb_num, Basis) #Quad C20
    
    rho_initial = np.zeros([12])
    rho_initial[:] =  0; rho_initial += 0.5
    rho_initial[[5, 7, 6]] = 10
    u00, u01, v_x, v_c, v_eff, rho_new = DFT_SCF(dual_faces, 0.001, rho_initial, orb_num, Basis1) #Quad C20
    
    fig = plt.figure()
    #fig.suptitle('Heat dissipation in C20')
    plt.subplot(221)
    #norm1 = matplotlib.colors.Normalize(vmin=20, vmax=100)
    for i in range(6):
        ax = plt.subplot(3,2,i+1)
        if i == 0:
            u = np.copy(u00)
            ax.set_title(r'(a) $\rho$')
        elif i == 1:
            u = np.copy(u01)
            ax.set_title(r'(b) $v_H$')
        elif i == 2:
            u = np.copy(u0)
            ax.set_title(r'(c) $\rho$')
        elif i == 3:
            u = np.copy(u1)
            ax.set_title(r'(d) $v_H$')
        elif i == 4:
            u = np.copy(u2)
            ax.set_title(r'(e) $\rho$')
        elif i == 5:
            u = np.copy(u3)
            ax.set_title(r'(f) $v_H$')
        data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
        if i > 1:
            A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
        else:
            A = barycentric_animation.contineous_rep(data_trans, u, Basis1, resolution, dual_faces)
        A1 = masked_array(A, A == -100)
    #A =np.copy(A1)

    #fig,ax = plt.subplots()
        if i == 1 or i == 3:
            cmap_alt = cm.bone
            pb = ax.imshow(A1, interpolation='nearest', cmap=cm.bone, vmax = 0.4)
        elif i == 5:
            cmap_alt = cm.bone
            pb = ax.imshow(A1, interpolation='nearest', cmap=cm.bone, vmax = 0.4)
        else:
            pb = ax.imshow(A1, interpolation='nearest', cmap=cm.inferno)
        
        #pb = ax.imshow(A1, interpolation='nearest', cmap=cmap_alt)
        
        #fig.colorbar(pb, cax=cbar, orientation='vertical')
        plt.colorbar(pb)
        plt.axis('off')
        #
    #fig.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #fig.colorbar(pb, cax=cbar_ax)
    plt.show()
    return
#import  FEM.BasisFunctions.BasisFunctions_Linear as Basis
import  FEM.BasisFunctions.BasisFunctions_Quadratic as Basis
import  FEM.BasisFunctions.BasisFunctions_Linear as Basis1
plot_pic_alt(Basis, dual_faces, arcs, arcpos, 400, Basis1)
#%%


def plot_pic_new(rho_full, con,Basis, dual_faces, arcs, arcpos, resolution):
    from Plotting import barycentric_animation;     from numpy.ma import masked_array
    from matplotlib import pyplot as plt;     from matplotlib import cm
    import matplotlib
    import matplotlib.gridspec as gridspec
    import matplotlib
    #fig = plt.figure()
    fig = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(3, 8)
    ax1 = plt.subplot(gs[0, 0:4])
    ax2 = plt.subplot(gs[0,4:])
    ax3 = plt.subplot(gs[1,0:4])
    ax4 = plt.subplot(gs[1,4:])
    ax5 = plt.subplot(gs[2,1:7])
    '''
    ax5 = plt.subplot(gs[2,0:4])
    ax6 = plt.subplot(gs[2,4:])
    
    ax7 = plt.subplot(gs[3,1:7])
    '''
    #ax6 = plt.subplot(gs[2,2:])
    #ax7 = plt.subplot(gs[3,1:3])
    #fig = gcf()
    gs.tight_layout(fig)
    ax_lst = [ax1,ax2,ax3,ax4,ax5]
    ind = [0, 1, 2, -1]
    for i in range(4):
        u = np.copy(rho_full[ind[i]])
        data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
        A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
        A1 = masked_array(A, A == -100)
        if i != 3:
            pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = np.max(rho_full[0]))
        elif i == 3:
            pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = np.min(u)-0.01, vmax = np.max(u)+0.01)
            
        '''
        elif i == 1:
            pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = np.max(rho_full[0]))
        else: 
            pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = np.max(rho_full[0]))
        '''
        #fig.colorbar(pb, cax=cbar, orientation='vertical')
        #plt.colorbar(pb)
        ax_lst[i].set_axis_off()
        if i == 0:
            ax_lst[i].set_title(r'(a) Initial density $\rho_0$')
        elif i == 1:
            ax_lst[i].set_title(r'(b) $\rho_1$')
        elif i == 2:
            ax_lst[i].set_title(r'(c) $\rho_2$')
        elif i == 3:
            #ax_lst[i].set_title(r'(c) $\rho_24$')
            #ax_lst[i].set_title(r'(d) Final density $\rho_{23}$')
            ax_lst[i].set_title(r'(d) Final density $\rho_{24}$')
        #plt.axis('off')
        fig.colorbar(pb, ax=ax_lst[i])
    x = np.linspace(1,len(con), len(con))
    ax5.plot(x, con, '.-')
    ax5.set_yscale('log');     ax5.set_xlabel('Iterations')
    ax5.set_ylabel(r'$L^2$-norm Difference')
    ax5.set_xlim([1,len(con)]);     ax5.set_ylim([0,max(con) +  1])
    ax5.set_title('(f) Convergence of SCF Loop')
    plt.show()
    return
plot_pic_new(rho_f, con,Basis, dual_faces, arcs, arcpos, 10)
#%%

def plot_pic_new_alt(orbs, energies, v, Basis, dual_faces, arcs, arcpos, resolution):
    from Plotting import barycentric_animation;     from numpy.ma import masked_array
    from matplotlib import pyplot as plt;     from matplotlib import cm
    import matplotlib
    import matplotlib.gridspec as gridspec
    import matplotlib
    #fig = plt.figure()
    fig = plt.figure(constrained_layout=True)
    gs = gridspec.GridSpec(3, 8)
    ax1 = plt.subplot(gs[0, 0:4])
    ax2 = plt.subplot(gs[0,4:])
    ax3 = plt.subplot(gs[1,0:4])
    ax4 = plt.subplot(gs[1,4:])
    #ax5 = plt.subplot(gs[2,1:7])
    ax5 = plt.subplot(gs[2,0:4])
    ax6 = plt.subplot(gs[2,4:])
    ax_lst = [ax1,ax2,ax3,ax4,ax5,ax6]
    #ind = [0, 1, 2, -1]
    ind = [0,1,3,10,18]
    #ind = [0,1,2,3,4,5]
    for i in range(-1,5):
        if i == -1:
            u = np.copy(-v)    
            cmapping = cm.bone
        else:
            u = np.copy(orbs[:,ind[i]])
            cmapping = cm.inferno
            
        data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
        A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
        A1 = masked_array(A, A == -100)
        #pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
        
        pb = ax_lst[i+1].imshow(A1, interpolation='nearest', cmap=cmapping)
        if i == -1:
            pb = ax_lst[0].imshow(A1, interpolation='nearest', cmap=cmapping, vmin = np.min(-v), vmax = np.max(-v)+1)

        ax_lst[i].set_axis_off()
        if i == -1:
            ax_lst[0].set_title(r'(a) Kohn-Sham Pontential')
            fig.colorbar(pb, ax=ax_lst[0])
        elif i == 0:
            ax_lst[i+1].set_title(r'(b) Orbital 1 Density')
        elif i == 1:
            ax_lst[i+1].set_title(r'(c) Orbital 2 Density')
        elif i == 2:
            ax_lst[i+1].set_title(r'(d) Orbital 4 Density')
        elif i == 3:
            ax_lst[i+1].set_title(r'(e) Orbital 11 Density')
        elif i == 4:
            ax_lst[i+1].set_title(r'(f) Orbital 19 Density')
        if i != -1:
            fig.colorbar(pb, ax=ax_lst[i+1])
        
    #fig.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #fig.colorbar(pb, cax=cbar_ax)
    print(energies[ind])
    plt.show()
    return
plot_pic_new_alt(orb, energies, v, Basis, dual_faces, arcs, arcpos, 400)

#%%
from Plotting import barycentric_animation;     from numpy.ma import masked_array
from matplotlib import pyplot as plt;     from matplotlib import cm
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib
resolution = 400
#fig = plt.figure()
fig = plt.figure(constrained_layout=True)
gs = gridspec.GridSpec(4, 3)
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[0, 1])
ax3 = plt.subplot(gs[0, 2])

ax4 = plt.subplot(gs[1, 0])
ax5 = plt.subplot(gs[1, 1])
ax6 = plt.subplot(gs[1, 2])

ax7 = plt.subplot(gs[2, 0])
ax8 = plt.subplot(gs[2, 1])
ax9 = plt.subplot(gs[2, 2])

ax10 = plt.subplot(gs[3, 0])
ax11 = plt.subplot(gs[3, 1])
ax12 = plt.subplot(gs[3, 2])


ax_lst = [ax1,ax2,ax3,ax4,ax5,ax6, ax7, ax8, ax9, ax10, ax11, ax12]

ax_lst[0].set_title(r'(a) Final Density $\rho_{final}$')
data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
A = barycentric_animation.contineous_rep(data_trans, rho, Basis, resolution, dual_faces)
A1 = masked_array(A, A == -100)
        #pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
pb = ax_lst[0].imshow(A1, interpolation='nearest', cmap=cm.inferno)
ax_lst[0].set_axis_off()
fig.colorbar(pb, ax=ax_lst[0])

subplottag = ['(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']
for i in range(1,12):
    u = np.copy(orb[:,i-1])
    
    data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
    A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
    A1 = masked_array(A, A == -100)
        #pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
    pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno)
    ax_lst[i].set_axis_off()

    #ax_lst[i].set_title('Orbital ' + str(i+1) + ' Density, E = ' + str(np.round(energies[indices[i+12]],3)) + ' Hartrees')
    ax_lst[i].set_title(subplottag[i-1] + ' Orbital ' + str(i) + ' Density, E = ' + str(np.round(energies[i-1],3)) + ' Hartrees')
    fig.colorbar(pb, ax=ax_lst[i])
    
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#fig.colorbar(pb, cax=cbar_ax)

#%%
colours = ['black', 'paleturquoise', 'lightskyblue', 'palegreen']
plt.figure()
x = np.linspace(1,30,30)
picked_rho = np.array([4, 5, 7, 11, 12, 13, 16, 21]) - 1
#for i in range(13):
for j in range(len(picked_rho)):
    i = picked_rho[j]
    l = len(full_conver[i])
    x = np.linspace(1,l,l)
    #if i > 9:
        #plt.plot(x, full_conver[i], label = str(i+1), color = colours[i-10],  linestyle='dashed', marker='.')
        #plt.scatter(l, full_conver[i][-1], color = colours[i-10])
    
    plt.plot(x, full_conver[i], linestyle='dashed', marker='.', label = str(i+1))
    plt.scatter(l, full_conver[i][-1])
    
plt.yscale('log')
plt.xlim([1,30]); plt.ylim([10**(-7), 10**2])
plt.legend(loc="upper right", ncol = 3, fontsize = 15)

plt.plot([15,15], [-5, 100], color='m')

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.text(15.75, 10.5, 'Density Mixing\nInitialized', fontsize=10,
        verticalalignment='top', color='m',  bbox=props)
plt.ylabel(r'$L^2$-norm Difference', fontsize = 15)
plt.xlabel('Iterations', fontsize = 20)
plt.xticks(np.linspace(1,30,30))
#ax.text(0.05, 0.95, 'Density Mixing Initiated', transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

#plt.gca().legend(loc='center left', bbox_to_anchor=(1, 0.5))
#%%
plt.figure()
x = np.linspace(0,30,31)
for i in range(25):
    l = len(full_conver[i])
    x = np.linspace(1,l,l)
    plt.plot(x, full_energy[i])
plt.yscale('symlog')


#%%
from Plotting import barycentric_animation;     from numpy.ma import masked_array
from matplotlib import pyplot as plt;     from matplotlib import cm
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib
resolution = 50
#fig = plt.figure()
fig = plt.figure(constrained_layout=True)
gs = gridspec.GridSpec(4, 3)
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[0, 1])
ax3 = plt.subplot(gs[0, 2])

ax4 = plt.subplot(gs[1, 0])
ax5 = plt.subplot(gs[1, 1])
ax6 = plt.subplot(gs[1, 2])

ax7 = plt.subplot(gs[2, 0])
ax8 = plt.subplot(gs[2, 1])
ax9 = plt.subplot(gs[2, 2])

ax10 = plt.subplot(gs[3, 0])
ax11 = plt.subplot(gs[3, 1])
ax12 = plt.subplot(gs[3, 2])


ax_lst = [ax1,ax2,ax3,ax4,ax5,ax6, ax7, ax8, ax9, ax10, ax11, ax12]

'''
ax_lst[0].set_title(r'(a) Final Density $\rho_{final}$')
data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
A = barycentric_animation.contineous_rep(data_trans, rho, Basis, resolution, dual_faces)
A1 = masked_array(A, A == -100)
        #pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
pb = ax_lst[0].imshow(A1, interpolation='nearest', cmap=cm.inferno)
ax_lst[0].set_axis_off()
'''
#index_rho = [1, 4, 9, 11, 12, 16, 21, 25] # 8 solutions?
index_rho = [1, 4, 9, 11, 12, 16, 21, 25] # 8 solutions?
subplottag = ['(a)','(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)']
for i in range(0, len(index_rho)):
    u = np.copy(full_rho[index_rho[i] - 1])
    
    data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
    A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
    A1 = masked_array(A, A == -100)
        #pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
    pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno)
    ax_lst[i].set_axis_off()

    #ax_lst[i].set_title('Orbital ' + str(i+1) + ' Density, E = ' + str(np.round(energies[indices[i+12]],3)) + ' Hartrees')
    ax_lst[i].set_title(subplottag[i] + ' ' +  str(index_rho[i]) + ' Orbital Density')
    fig.colorbar(pb, ax=ax_lst[i])
    
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
#fig.colorbar(pb, cax=cbar_ax)
#%%
from Plotting import barycentric_animation;     from numpy.ma import masked_array
from matplotlib import pyplot as plt;     from matplotlib import cm
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib
resolution = 500
fig = plt.figure()
u = np.copy(rho)
    
data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
A1 = masked_array(A, A == -100)
#pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
pb = plt.imshow(A1, interpolation='nearest',cmap=cm.inferno, origin='lower')
fig.colorbar(pb, orientation='vertical')
plt.axis('off')
plt.title(r'Final Density $\rho_{final}$', fontsize = 15)

#%%
resolution = 50
orby = 30
fig = plt.figure()
u = np.copy(orb[:,orby-1])
data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
A1 = masked_array(A, A == -100)
#pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
pb = plt.imshow(A1, interpolation='nearest',cmap=cm.inferno, origin='lower')
fig.colorbar(pb, orientation='vertical')
plt.axis('off')
plt.title('Orbital ' + str(orby) + ' Density, E = ' + str(np.round(energies[orby-1],3)) + ' Hartrees', fontsize = 15)

print('Sim # finished')

orby = 31
fig = plt.figure()
u = np.copy(orb[:,orby-1])
data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
A1 = masked_array(A, A == -100)
#pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
pb = plt.imshow(A1, interpolation='nearest',cmap=cm.inferno, origin='lower')
fig.colorbar(pb, orientation='vertical')
plt.axis('off')
plt.title('Orbital ' + str(orby) + ' Density, E = ' + str(np.round(energies[orby-1],3)) + ' Hartrees', fontsize = 15)
                       
print('Sim # finished')                                         
                                                                
orby = 52
fig = plt.figure()
u = np.copy(orb[:,orby-1])
data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
A1 = masked_array(A, A == -100)
#pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
pb = plt.imshow(A1, interpolation='nearest',cmap=cm.inferno, origin='lower')
fig.colorbar(pb, orientation='vertical')
plt.axis('off')
plt.title('Orbital ' + str(orby) + ' Density, E = ' + str(np.round(energies[orby-1],3)) + ' Hartrees', fontsize = 15)

print('Sim # finished')

orby = 57
fig = plt.figure()
u = np.copy(orb[:,orby-1])
data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
A1 = masked_array(A, A == -100)
#pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
pb = plt.imshow(A1, interpolation='nearest',cmap=cm.inferno, origin='lower')
fig.colorbar(pb, orientation='vertical')
plt.axis('off')
plt.title('Orbital ' + str(orby) + ' Density, E = ' + str(np.round(energies[orby-1],3)) + ' Hartrees', fontsize = 15)

print('Sim # finished')

orby = 65
fig = plt.figure()
u = np.copy(orb[:,orby-1])
data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
A1 = masked_array(A, A == -100)
#pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
pb = plt.imshow(A1, interpolation='nearest',cmap=cm.inferno, origin='lower')
fig.colorbar(pb, orientation='vertical')
plt.axis('off')
plt.title('Orbital ' + str(orby) + ' Density, E = ' + str(np.round(energies[orby-1],3)) + ' Hartrees', fontsize = 15)
#%%       
from Plotting import barycentric_animation;     from numpy.ma import masked_array
from matplotlib import pyplot as plt;     from matplotlib import cm
import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib
import time
resolution = 400
fig = plt.figure(constrained_layout=True)
gs = gridspec.GridSpec(2, 2)
ax1 = plt.subplot(gs[0, 0]); ax2 = plt.subplot(gs[0, 1])
ax3 = plt.subplot(gs[1, 0]); ax4 = plt.subplot(gs[1, 1])

ax_lst = [ax1,ax2,ax3,ax4]
#orb_list = [35, 55, 57, 65]
orb_list = [31, 51, 55, 70]
for i in range(4):
    t0 = time.time()
    
    u = np.copy(orb[:,orb_list[i]-1])
    data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
    A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
    A1 = masked_array(A, A == -100)
        #pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 0.26)
    pb = ax_lst[i].imshow(A1, interpolation='nearest', cmap=cm.inferno, origin='lower')
    ax_lst[i].set_axis_off()
    
    #ax_lst[i].set_title('Orbital ' + str(i+1) + ' Density, E = ' + str(np.round(energies[indices[i+12]],3)) + ' Hartrees')
    ax_lst[i].set_title('Orbital ' + str(orb_list[i]) + ' Density, E = ' + str(np.round(energies[orb_list[i]-1],3)) + ' Hartrees', fontsize = 15)
    fig.colorbar(pb, ax=ax_lst[i])
    print('Plot ' + str(i) + ' computed')
    
    t1 = time.time()
    total = t1-t0
    print(total/60)
"Approximately 6 hours of run-time for these 4 orbitals at a resolution of 400 pixels across the x-axis"

#%%
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.special import sph_harm

phi = np.linspace(0, np.pi, 100)
theta = np.linspace(0, 2*np.pi, 100)
phi, theta = np.meshgrid(phi, theta)

# The Cartesian coordinates of the unit sphere
x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)

m, l = 1, 1

Y__1_1 = 1/2 * np.sqrt(3/(2*np.pi)) * np.e**(-1j*theta) * np.sin(phi)
Y_0_1 = 1/2 * np.sqrt(3/(np.pi)) * np.cos(phi)
Y_1_1 = -1/2 * np.sqrt(3/(2*np.pi))* np.e**(1j*(theta + np.pi/2)) * np.sin(phi)

# Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
#fcolors = (sph_harm(1, l, theta, phi).real)
#fcolors = (sph_harm(-1, l, theta, phi).real)**2 + (sph_harm(1, l, theta, phi).real)**2 #+ (sph_harm(0, l, theta, phi).real)**2
fcolors = (Y_0_1.real)**2 
#fcolors = (Y__1_1 * np.conj(Y_1_1)).real
        
          
#fcolors = sph_harm(-1, l, theta, phi).real + sph_harm(0, l, theta, phi).real + sph_harm(1, l, theta, phi).real
fmax, fmin = fcolors.max(), fcolors.min()
fcolors = (fcolors - fmin)/(fmax - fmin)

# Set the aspect ratio to 1 so our sphere looks spherical
fig = plt.figure(figsize=plt.figaspect(1.))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors))
# Turn off the axis planes
ax.set_axis_off()
plt.show()



#%%
'''
    m =  0
    for i in range(len(triangulation)):
        m1 = []
        m1.append(np.max( np.abs(rho_old[all_nodes[i]] - np.roll(rho_old[all_nodes[i]], 1))))
        m1.append(np.max( np.abs(rho_old[all_nodes[i]] - np.roll(rho_old[all_nodes[i]], 2))))
        m1.append(np.max( np.abs(rho_old[all_nodes[i]] - np.roll(rho_old[all_nodes[i]], 3))))
        m1.append(np.max( np.abs(rho_old[all_nodes[i]] - np.roll(rho_old[all_nodes[i]], 4))))
        m1.append(np.max( np.abs(rho_old[all_nodes[i]] - np.roll(rho_old[all_nodes[i]], 5))))
        if max(m1) > m:
            m = max(m1)
    m_reference = m
    #print('m reference is ' + str(m))
    if m_reference == 0:
        m_reference = 3
        '''
