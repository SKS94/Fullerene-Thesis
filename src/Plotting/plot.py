#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 15:36:59 2020

@author: simon
"""

def plot_2D(rho, basis, resolution, title, dual_faces, arcs, arcpos):
    import numpy as np    
    from matplotlib import pyplot as plt
    #grid, arcs, arcpos = data()

    #dual_faces = np.array([[6,7,9],[6,9,8],[5,7,6],[7,11,9],[4,7,5],[3,6,8],[2,8,10],[2,3,8],[9,11,10],[8,9,10],[4,11,7],[3,5,6],[1,10,11],[1,2,10],[1,11,4],[0,4,5],[0,1,4],[0,2,1],[0,3,2],[0,5,3]])
    fig,ax = plt.subplots()
    
    trans_mat = np.array([[1, np.cos(np.pi/3)],[0, np.sin(np.pi/3)]])
    
    from numpy.ma import masked_array
    import Plotting.barycentric_animation as barycentric_animation
    data_trans = barycentric_animation.transformedFromEisenstein(dual_faces, arcs, arcpos)

    A = barycentric_animation.contineous_rep(data_trans, rho, basis, resolution, dual_faces)
    A1 = masked_array(A, A == -100)
    #A =np.copy(A1)
    from matplotlib import cm
    #fig,ax = plt.subplots()

    #pb = ax.imshow(A1, interpolation='nearest',cmap=cm.inferno,  vmin=-0.05, vmax=0.6296296296296299)
#    pb = ax.imshow(A1, interpolation='nearest',cmap=cm.inferno, vmin = 0, vmax = 0.04)
    pb = ax.imshow(A1, interpolation='nearest',cmap=cm.RdBu, origin='upper', vmin = -0.35, vmax = 0.35)
    fig.colorbar(pb, orientation='vertical')
    
    plt.title(title)
    plt.show()
    return pb


def data():
    points = 12
    grid = [[-1, 10, 10, 10, 10, 10],[2, 8, 9, 11, 1, 2],[3, 6, 7, 4, 0, 3],[5, 5, 5, 5, 5, -1]];

    arcs = [[2, 3], [3, 5], [6, 5], [7, 5], [4, 5], [0, 5], [5, 6], [5, 7], [5, 4], [5, 0], [5, 3], [8, 10], [9, 10], [11, 10], [1, 10], 
            [2, 10], [3, 2], [10, 2], [10, 8], [10, 9], [10, 11], [10, 1], [2, 8], [8, 9], [9, 11], [11, 1], [1, 2], [3, 6], [6, 7], [7, 4],
            [4, 0], [0, 3], [8, 6], [9, 7], [11, 4], [1, 0], [3, 8], [6, 9], [7, 11], [4, 1], [0, 2], [8, 2], [9, 8], [11, 9], [1, 11], [2, 1],
            [6, 3], [7, 6], [4, 7], [0, 4], [3, 0], [6, 8], [7, 9], [4, 11], [0, 1], [8, 3], [9, 6], [11, 7], [1, 4], [2, 0]];
    arcpos = [[[2, 1], [3, 1]], [[3, 1], [4, 1]], [[3, 2], [4, 2]], [[3, 3], [4, 3]], [[3, 4], [4, 4]], [[3, 5], [4, 5]], [[4, 1], [3, 2]],
              [[4, 2], [3, 3]], [[4, 3], [3, 4]], [[4, 4], [3, 5]], [[4, 5], [3, 6]], [[2, 2], [1, 2]], [[2, 3], [1, 3]], [[2, 4], [1, 4]],
              [[2, 5], [1, 5]], [[2, 6], [1, 6]], [[3, 6], [2, 6]], [[1, 2], [2, 1]], [[1, 3], [2, 2]], [[1, 4], [2, 3]], [[1, 5], [2, 4]],
              [[1, 6], [2, 5]], [[2, 1], [2, 2]], [[2, 2], [2, 3]], [[2, 3], [2, 4]], [[2, 4], [2, 5]], [[2, 5], [2, 6]], [[3, 1], [3, 2]],
              [[3, 2], [3, 3]], [[3, 3], [3, 4]], [[3, 4], [3, 5]], [[3, 5], [3, 6]], [[2, 2], [3, 2]], [[2, 3], [3, 3]], [[2, 4], [3, 4]],
              [[2, 5], [3, 5]], [[3, 1], [2, 2]], [[3, 2], [2, 3]], [[3, 3], [2, 4]], [[3, 4], [2, 5]], [[3, 5], [2, 6]], [[2, 2], [2, 1]],
              [[2, 3], [2, 2]], [[2, 4], [2, 3]], [[2, 5], [2, 4]], [[2, 6], [2, 5]], [[3, 2], [3, 1]], [[3, 3], [3, 2]], [[3, 4], [3, 3]],
              [[3, 5], [3, 4]], [[3, 6], [3, 5]], [[3, 2], [2, 2]], [[3, 3], [2, 3]], [[3, 4], [2, 4]], [[3, 5], [2, 5]], [[2, 2], [3, 1]],
              [[2, 3], [3, 2]], [[2, 4], [3, 3]], [[2, 5], [3, 4]], [[2, 6], [3, 5]]];
              

    return grid, arcs, arcpos
    