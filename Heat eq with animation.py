#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:20:59 2020
@author: simon


Solving the Heat Equation using FEM on the C20 Fullerene.
Forward Euler is used for time evolution.


! With added animation to visualize the solution.
"""
#%%
from matplotlib import pyplot as plt
#fig = plt.figure(1)
fig, ax = plt.subplots()

ax.set_xlim(0,6.5); ax.set_ylim(-0.5,3)
ax.set_xticks([]);  ax.set_yticks([])

import numpy as np
Einstein_grid = np.zeros([4,6], dtype = int)
Einstein_grid[0,:] = 10; Einstein_grid[0,0] = -1
Einstein_grid[-1,:] = 5; Einstein_grid[-1,-1] = -1
Einstein_grid[1,0] = 2; Einstein_grid[1,1] = 8; Einstein_grid[1,2] = 9; Einstein_grid[1,3] = 11; Einstein_grid[1,4] = 1; Einstein_grid[1,5] = 2;
Einstein_grid[2,0] = 3; Einstein_grid[2,1] = 6; Einstein_grid[2,2] = 7; Einstein_grid[2,3] = 4; Einstein_grid[2,4] = 0; Einstein_grid[2,5] = 3;
print(Einstein_grid)

coordinate_grid = np.zeros([])
trans_mat = np.array([[1, np.cos(np.pi/3)],[0, np.sin(np.pi/3)]])

twoD_coor = {}
for i in range(Einstein_grid.shape[0]):
    for j in range(Einstein_grid.shape[1]):
        if Einstein_grid[i,j] != -1:
            
            coor = np.array([[j],[i]])            
            equil_coor = np.dot(trans_mat, coor).T[0]
            scatter = ax.scatter([equil_coor[0], equil_coor[0]], [equil_coor[1], equil_coor[1]], c = [0,0], cmap=plt.cm.plasma , s = 150)
            if str(Einstein_grid[i,j]) in twoD_coor:
                twoD_coor[str(Einstein_grid[i,j])] = np.vstack((twoD_coor[str(Einstein_grid[i,j])],   equil_coor))
            else:
                twoD_coor[str(Einstein_grid[i,j])] = equil_coor

#%%
import FEM.FEM_Assembly as FEM_Assembly
from FEM.TriangularElements.shape_construction import FE_construction
import  FEM.BasisFunctions.BasisFunctions_Quadratic as Basis

import os.path

import numpy as np
dual_faces = np.array([[6,7,9],[6,9,8],[5,7,6],[7,11,9],[4,7,5],[3,6,8],[2,8,10],[2,3,8],[9,11,10],[8,9,10],[4,11,7],[3,5,6],[1,10,11],[1,2,10],[1,11,4],[0,4,5],[0,1,4],[0,2,1],[0,3,2],[0,5,3]])

folder = 'Coordinate_Quadrature_data'
coordinates_7 = np.load(os.path.join(folder, 'coordinates_7.npy'))
weights_7 = np.load(os.path.join(folder, 'weights_7.npy'))

#%%
'''
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

'''
#W = Stiffness matrix,       M = Mass matrix,       b = Load vector
'''
#W, b, M = FEM_Assembly.Assemble(dual_faces, Element, BasisFunctions_Quadratic, coordinates_7, weights_7)
W, b, M = FEM_Assembly.Assemble(dual_faces, Element, BasisFunctions_Linear, coordinates_7, weights_7)

u_old = np.zeros([len(b)]);
delta_time = 0.001;     steps = 10000


for i in range(steps):        
    lin_sys_vec = np.dot(M, u_old) - delta_time * np.dot(W, u_old) +  b
    u_new = np.linalg.solve(M, lin_sys_vec)
    
    u_old = np.copy(u_new)  #Setting the heat source term to zero after the first step. 
    if i == 0:
        b[:] = 0

    for j in range(12):
        value_sum = np.sum(u_new)
        value = u_new[j] 
        t = value/value_sum
        for k in range(twoD_coor[str(j)].shape[0]):
            try:
                a = twoD_coor[str(j)].shape[1]
                coor = twoD_coor[str(j)][k]
                plt.scatter(coor[0], coor[1], c = t)
            except IndexError:
                coor = twoD_coor[str(j)]
                plt.scatter(coor[0], coor[1], c = t)
                
    if i//100 == i / 100:
        plt.show()
        print(i)
'''
#%%
'''
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import cm

x = []; y = []        
fig, ax = plt.subplots()


ax.set_xlim(0,105)
ax.set_ylim(0,12)
line, = ax.plot(0,0)

def ani(i):
    x.append(i * 10)
    y.append(i)
    
    line.set_xdata(x)
    line.set_ydata(y)
    return line


animation = FuncAnimation(fig, func = ani, frames = np.arange(0, 100, 0.01), interval = 200)
plt.show()
'''
#%%
f = np.zeros([42]) + 0.00001
f[6] = 1
M, W, b = FEM_Assembly.Assemble(dual_faces, Basis, coordinates_7, weights_7, 'Heat', f, 'dual')
M = M.toarray();
W = W.toarray();
#u_old = np.zeros([len(b)]);
u_old = np.zeros([42])


delta_time = 0.01;     steps = 100

color_val = [np.sum(b)/3, np.sum(b)/3, np.sum(b)/3]

x = [twoD_coor['6'][0], twoD_coor['7'][0], twoD_coor['9'][0]]
y = [twoD_coor['6'][1], twoD_coor['7'][1], twoD_coor['9'][1]]

ax.scatter(x, y,  c=color_val, cmap=plt.cm.plasma, s = 150, vmin=0, vmax=1/3)
ax.set_title('Heat dissipation C20')

lin_sys_vec = np.dot(M, u_old) - delta_time * np.dot(W, u_old) +  b
u_new = np.linalg.solve(M, lin_sys_vec)
    
u_old = np.copy(u_new)  #Setting the heat source term to zero after the first step. 

#b[:] = 0
#%%
import time
start = time.time()

import matplotlib.animation as animation
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

from matplotlib.animation import FuncAnimation

def heat_ani(i):
    global u_old; global u_new
    
    scatter_points_x = [];  scatter_points_y = []; color_val = []
    lin_sys_vec = np.dot(M, u_old) - delta_time * np.dot(W, u_old) +  b
    u_new = np.linalg.solve(M, lin_sys_vec)
    
    u_old = np.copy(u_new)  #Setting the heat source term to zero after the first step. 

    for j in range(12):
        value_sum = np.sum(u_new)
        value = u_new[j] 
        t = value/value_sum
            
        for k in range(twoD_coor[str(j)].shape[0]):
            try:
                a = twoD_coor[str(j)].shape[1]
                coor = twoD_coor[str(j)][k]

                scatter_points_x.append(coor[0]);   scatter_points_y.append(coor[1])
                color_val.append(t)
            except IndexError:
                coor = twoD_coor[str(j)]
                scatter_points_x.append(coor[0]);   scatter_points_y.append(coor[1])
                color_val.append(t)
                
    #print(np.sum(u_new))
    ax.scatter(scatter_points_x, scatter_points_y, c=color_val, cmap=plt.cm.plasma, s = 150, vmin=-0.1, vmax=1/3)
    #plt.pause(0.01)
    return scatter
                
animation1 = FuncAnimation(fig, func = heat_ani, frames = np.arange(0, 100, 1), interval = 1, repeat = False, save_count=200)
#%%
#FFwriter = animation.FFMpegWriter(fps=30, codec="libx264")  
#writer = animation.writers['avconv'](fps=30)
#writer = animation.writers['ffmpeg']
#%%
#animation1.save('basic_animation1.mp4', writer = writer )
animation1.save('heat5.gif', dpi=80, writer='imagemagick', fps=30)

plt.show()
end = time.time()
print('Time for animation:' +  str(end - start))
#%%
#fig, ax = plt.subplots()
#heat_ani(2)



