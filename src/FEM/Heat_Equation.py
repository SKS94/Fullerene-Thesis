#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:20:59 2020

@author: simon


Solving the Heat Equation using FEM on the C20 Fullerene.
Forward Euler is used for time evolution.
"""

'''
import FEM_Assembly
import BasisFunctions_Quadratic
import BasisFunctions_Linear
from shape_construction import FE_construction
import os.path

import numpy as np
#dual_faces = np.array([[6,7,9],[6,9,8],[5,7,6],[7,11,9],[4,7,5],[3,6,8],[2,8,10],[2,3,8],[9,11,10],[8,9,10],[4,11,7],[3,5,6],[1,10,11],[1,2,10],[1,11,4],[0,4,5],[0,1,4],[0,2,1],[0,3,2],[0,5,3]])
triangles_shape = [60,3];
dual_faces = np.array([[[0,16,15],[0,15,14],[0,14,13],[0,13,12],[0,12,16],[1,17,22],[1,22,21],[1,21,16],[1,16,12],[1,12,17],[2,13,18],[2,18,23],[2,23,17],[2,17,12],[2,12,13],[3,14,19],[3,19,24],[3,24,18],[3,18,13],[3,13,14],[4,15,20],[4,20,25],[4,25,19],[4,19,14],[4,14,15],[5,16,21],[5,21,26],[5,26,20],[5,20,15],[5,15,16],[6,23,28],[6,28,27],[6,27,22],[6,22,17],[6,17,23],[7,24,29],[7,29,28],[7,28,23],[7,23,18],[7,18,24],[8,25,30],[8,30,29],[8,29,24],[8,24,19],[8,19,25],[9,26,31],[9,31,30],[9,30,25],[9,25,20],[9,20,26],[10,22,27],[10,27,31],[10,31,26],[10,26,21],[10,21,22],[11,31,27],[11,27,28],[11,28,29],[11,29,30],[11,30,31]]]).reshape(triangles_shape);
Element = FE_construction('triangle_equil')
folder = 'Coordinate_Quadrature_data'
coordinates_7 = np.load(os.path.join(folder, 'coordinates_7.npy'))
weights_7 = np.load(os.path.join(folder, 'weights_7.npy'))
'''
#%%
'''
W = Stiffness matrix,       M = Mass matrix,       b = Load vector
'''


'''
W, b, M = FEM_Assembly.Assemble(dual_faces, Element, BasisFunctions_Quadratic, coordinates_7, weights_7)
#W, b, M = FEM_Assembly.Assemble(dual_faces, Element, BasisFunctions_Linear, coordinates_7, weights_7)

u_old = np.zeros([len(b)]);
u_old = np.copy(b)
b[:] = 0

delta_time = 0.05
steps = 200

for i in range(steps):        
    #lin_sys_vec = np.dot(M, u_old) - delta_time * np.dot(W, u_old) +  b
    lin_sys_vec = np.dot(M, u_old) - delta_time * np.dot(W, u_old) +  b
    u_new = np.linalg.solve(M, lin_sys_vec)
    
    
    div = np.abs(np.abs(u_old[0]) - np.abs(u_new[0]))
    if div < 0.000000001:
        break
    
    u_old = np.copy(u_new)  #Setting the heat source term to zero after the first step. 
    if i == 0:
        b[:] = 0
        
    if i < 30:
        print(np.max(u_new))
'''
#%%    

