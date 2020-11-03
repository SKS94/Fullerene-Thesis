#%%
import numpy as np
import os.path
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.close('all')
from matplotlib.animation import FuncAnimation
from Plotting import barycentric_animation

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



import FEM.FEM_Assembly as FEM_Assembly
from FEM.TriangularElements.shape_construction import FE_construction
#import  FEM.BasisFunctions.BasisFunctions_Quadratic as Basis
import  FEM.BasisFunctions.BasisFunctions_Linear as Basis

points = Basis.GlobalDOF(dual_faces)
import os.path
folder = 'Coordinate_Quadrature_data'
coordinates_7 = np.load(os.path.join(folder, 'coordinates_7.npy'))
weights_7 = np.load(os.path.join(folder, 'weights_7.npy'))
resolution = 10
#%%
for i in range(len(arcpos)):
    x = [arcpos[i][0][0], arcpos[i][1][0]]
    y = [arcpos[i][0][1], arcpos[i][1][1]]
    plt.plot(y,x, 'b')
plt.show()

#%%
def colormap_test(M, W, u_test, delta_time, steps):
    maximum = np.max(u_test);   minimum = np.min(u_test)
    for i in range(steps):        
        lin_sys_vec = np.dot(M, u_test) - delta_time * np.dot(W, u_test) + b*delta_time
        u_test_new = np.linalg.lstsq(M, lin_sys_vec)[0]
        u_test = np.copy(u_test_new)  #Setting the heat source term to zero after the first step. 
        if np.max(u_test) > maximum:
            maximum = np.max(u_test)
        if np.min(u_test) < minimum:
            minimum = np.min(u_test)
    return maximum, minimum

#%%
'''f = np.zeros([42]) + 0.00001
f = np.ones([42])*0.5
f[6] = 10
f[7] = 10
f[9] = 10
f[12] = 10;  f[13] = 10;  f[14] = 10'''
'''f = np.zeros([42]) + 0.00001
f = np.ones([42])*0.5
f[6] = 10'''
delta_time = 0.015; steps = 80
#delta_time = 0.01;   steps = 350
trans_mat = np.array([[1, np.cos(np.pi/3)],[0, np.sin(np.pi/3)]])

#p = np.random.rand(42)
#p[:] *= 100

#p = np.zeros([42])
#p[:] += 50
#p[5] = 100

p = np.array([0., 0., 0., 0., 0., 0., 1., 1., 0., 1., 0., 0., 1., 1., 1., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0.]) * 100

M, W, b = FEM_Assembly.Assemble(dual_faces, Basis, coordinates_7, weights_7, 'Heat', p, 'dual')
#M, W, b = FEM_Assembly.Assemble(dual_faces, Basis, coordinates_7, weights_7, 'Heat', p, 'dual')
M = M.toarray();    W = W.toarray();

#u_old = np.zeros([points]);    u_old = np.copy(b);
b=0
u_old = np.copy(p)
#u_old[10] = 5
#u_old = np.zeros([42]); u_old[10] = 20
maximum, minimum = colormap_test(M, W, u_old, delta_time, steps)


#%%
import string;    import random
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

name = id_generator()
file1 = open('HeatGifs/' + str(Fullerene) + '_' + str(name) + ".txt","a") 
#maximum = 5; minimum = 0
text = str('Number of steps ') + str(steps) + ' using a delta time of ' + str(delta_time) +str(' at a resolution of ') + str(resolution) + \
           '\n with an initial u of ' + str(u_old)
file1.write(text)

#fig, ax = plt.subplots()

fig,[ax,cax] = plt.subplots(1,2, gridspec_kw={"width_ratios":[50,1]})
# Set the colormap and norm to correspond to the data for which5
# the colorbar will be used.
from matplotlib import cm
import matplotlib
cmap = matplotlib.cm.inferno
norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)
#norm = matplotlib.colors.LogNorm()
cb1 = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')

ax.set_title('Heat dissipation in C20 0 iter')
plt.axis('off')
#%%
def ani(i):
    global u_new; global u_old; global b
    if i == 0:
        plot_u(i, u_old)
        return
    if i < 20:
        plot_u(i, u_old)
        return
#    lin_sys_vec = np.dot(M, u_old) - delta_time * (np.dot(np.dot(np.linalg.inv(M),W), u_old)) +  delta_time * b
    lin_sys_vec = np.dot(M, u_old) - delta_time * (np.dot(W, u_old)) +  delta_time * b
    u_new = np.linalg.solve(M, lin_sys_vec)
    u_old = np.copy(u_new)  #Setting the heat source term to zero after the first step. 
    plot_u(i, u_old)
    if i == 50:
        b = 0 
    
    return
     
rho_0 = np.zeros([points]); rho_1 = np.zeros([points]); rho_2 = np.zeros([points])
rho_3 = np.zeros([points])
def plot_u(i, u):
    import time
    start = time.time()
    
    from Plotting import barycentric_animation
    global points; global trans_mat
    global maximum, minimum, rho_0, rho_1, rho_2, rho_3;
    global dual_faces, arcs, arcpos, grid, resolution, u_old;
    print(i)
    
    if i == 0:
        rho_0 = np.copy(u)
    elif i == 25:
        rho_1 = np.copy(u)
    elif i == 40:
        rho_2 = np.copy(u)
    elif i == 70:
        rho_3 = np.copy(u)
#    basis = BasisFunctions_Quadratic

    from numpy.ma import masked_array

    data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)

    A = barycentric_animation.contineous_rep(data_trans, u, Basis, resolution, dual_faces)
    A1 = masked_array(A, A == -100)
    #A =np.copy(A1)
    from matplotlib import cm
    #fig,ax = plt.subplots()

    pb = ax.imshow(A1, interpolation='nearest', cmap=cm.inferno,  vmin=minimum, vmax=maximum)
    #if i == 0:
        #cbar = ax.colorbar(pb)
    #pb = ax.imshow(A1, interpolation='nearest',cmap=cm.inferno,  vmin=-0.2, vmax=0.6296296296296299)
    
    if i < 20:
        ax.set_title('Heat dissipation in C20 at iteration: ' + str(0))
    else:
        ax.set_title('Heat dissipation in C20 at iteration: ' + str(i - 20))
        
    end = time.time()
    print('Simulation and animation took: ' + str((end - start))+ ' s')
    return pb
ani(11)
#a = plot_u(1, u_old)
#plt.show()
#%%
plt.show()

import time
start = time.time()
#animation1 = FuncAnimation(fig, func = ani, frames = np.arange(0, 350, 1), interval = 0, repeat = False, save_count=350)
animation1 = FuncAnimation(fig, func = ani, frames = np.arange(0, 100, 1), interval = 0, repeat = False, save_count=200)
#fig.colorbar(pb, cax=cax, orientation='vertical')
#plt.axis('off')
#ax.set_xticks([]); ax.set_yticks([])
plt.show()
#%%
animation1.save('HeatGifs/' + Fullerene + '_' + name + '.gif', writer='imagemagick', dpi=80, fps=10)

end = time.time()

print('Simulation and animation took: ' + str((end - start)/60)+ ' min.')
file1.close()
#%%
def plot_pic(u0, u1, u2, u3, Basis, dual_faces, arcs, arcpos, resolution, a, b, c, d):
    from numpy.ma import masked_array
    from matplotlib import pyplot as plt
    from matplotlib import cm
    import matplotlib
    fig = plt.figure(3)
    fig.suptitle('Heat dissipation in C20')
    plt.subplot(221)
    norm1 = matplotlib.colors.Normalize(vmin=20, vmax=100)
    for i in range(4):
        ax = plt.subplot(2,2,i+1)
        if i == 0:
            u = np.copy(u0)
            num = a - 20
            ax.set_title('(a) Iteration: ' + str(num))
        elif i == 1:
            u = np.copy(u1)
            num = b - 20
            ax.set_title('(b) Iteration: ' + str(num))
        elif i == 2:
            u = np.copy(u2)
            num = c - 20
            ax.set_title('(c) Iteration: ' + str(num))
        elif i == 3:
            u = np.copy(u3)
            num = d - 20
            ax.set_title('(d) Iteration: ' + str(num))
        
        
        data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)
        A = barycentric_animation.contineous_rep(data_trans, u, Basis, 400, dual_faces)
        A1 = masked_array(A, A == -100)
    #A =np.copy(A1)
    
    #fig,ax = plt.subplots()
    
        pb = ax.imshow(A1, interpolation='nearest', cmap=cm.inferno, vmin = 0, vmax = 100)
        plt.axis('off')
        #
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(pb, cax=cbar_ax)

    plt.show()
    return


#%%
delta_time = 0.01
u_old = rho_0
for i in range(500):
    lin_sys_vec = np.dot(M, u_old) - delta_time * (np.dot(W, u_old)) +  delta_time * b
    u_new = np.linalg.solve(M, lin_sys_vec)
    u_old = np.copy(u_new)  #Setting the heat source term to zero after the first step. 
    if i == 30:
        rho_1 = np.copy(u_old)
    elif i == 100:
        rho_2 = np.copy(u_old)
rho_3 = np.copy(u_new)
#%%
plot_pic(rho_0, rho_1, rho_2, rho_3, Basis, dual_faces, arcs, arcpos, resolution, 20, 50, 120, 520)
