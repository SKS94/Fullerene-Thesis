#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 14:04:34 2020

@author: simon

Torus generate
"""
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import array
plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.set_xlim([0,2])
#%% C168
coor = np.zeros([168,3]);   count = 7

fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')

coor[:7,:] = array([[0.664, -4.344, 0.666], [1.569, -4.860, 1.703], [1.334, -6.284, 1.890], [0, -6.875, 1.625], [0.733, -7.364, -1.537],
              [1.167, -7.960, -0.255], [0, -7.859, 0.617]])
ax.scatter(coor[:7,0], coor[:7,1], coor[:7,2], s=20, c='g')
reflect_yz = np.identity(3);    reflect_yz[0,0] = -1
num = [0,1,2,4,5]
for ind in num:
    vec = (coor[ind,:]).T 

    new = np.dot(reflect_yz, vec)
    coor[count,:] = new
    ax.scatter(new[0], new[1], new[2], s=20, c='g')
    count += 1

#%%
from numpy import cos, sin
u = np.cos(2*np.pi/7)
v = np.sin(2*np.pi/7)
w = 0 ; a=u; b=v; c = 0
angle = np.pi
for i in range(count):    
    vec =  coor[i,:]
    x = vec[0]; y =  vec[1];  z = vec[2]
#    new = np.dot(rot, vec)
    #coor[count,:] = new
    xn = (a*v**2-u*(b*v-u*x-v*y))*(1-cos(angle))+x*cos(angle)
    yn = (b*u**2-v*(a*u-u*x-v*y))*(1-cos(angle))+y*cos(angle)
    zn = (c*(u**2+v**2)) + z*cos(angle) +(-b*u + a*v - v*x+ u*y)*sin(angle)
    ax.scatter(xn, yn, zn, s=20, c='g')
    coor[count,0] =  xn;    coor[count,1] =  yn;    coor[count,2] =  zn
    count += 1

#%%
rot_cn_z = array([[np.cos(2*np.pi/5), -np.sin(2*np.pi/7), 0],
            [np.sin(2*np.pi/5), np.cos(2*np.pi/7),  0],
            [0,0,1]])

for i in range(0, count):
    vec = (coor[i,:]).T 
    for j in range(1,7):
        rot_cn_z = array([[np.cos(j*2*np.pi/7), -np.sin(j*2*np.pi/7), 0],
            [np.sin(j*2*np.pi/7), np.cos(j*2*np.pi/7),  0],
            [0,0,1]])
        new = np.dot(rot_cn_z, vec)
        coor[count,:] = new
        ax.scatter(new[0], new[1], new[2], s=20, c='g')
        count+=1
#%%
for i in range(168):
    ax.scatter(coor[i,0], coor[i,1], coor[i,2], s=20, c='g')
    diff = (coor - coor[i,:])
    norm = np.linalg.norm(diff,axis=1)
    idx = np.argpartition(norm, 4)
    tick = 0
    while True:
        if norm[idx[tick]] > 0.1 :
            x = [coor[idx[tick],0], coor[i,0]]
            y = [coor[idx[tick],1], coor[i,1]]
            z = [coor[idx[tick],2], coor[i,2]]
            ax.plot(x,y,z,'y')
        tick+=1
        if tick == 4:
            break
#%% Norm based sorting of polygons
for i in range(168):
    diff = (coor - coor[i,:])
    norm = np.linalg.norm(diff,axis=1)
    idx = np.argpartition(norm, 4)
    tick = 0
    while True:
        if 1.3 < norm[idx[tick]] < 1.4 :
            x = [coor[idx[tick],0], coor[i,0]]
            y = [coor[idx[tick],1], coor[i,1]]
            z = [coor[idx[tick],2], coor[i,2]]
            ax.plot(x,y,z,c='b')
        elif 1.4 < norm[idx[tick]] < 1.5:
            x = [coor[idx[tick],0], coor[i,0]]
            y = [coor[idx[tick],1], coor[i,1]]
            z = [coor[idx[tick],2], coor[i,2]]
            ax.plot(x,y,z,c='r')
        tick+=1
        if tick == 4:
            break
#%%

def neighbour(coor, point):
    new_coor = np.zeros([3,3])
    diff = (coor - point)
    norm = np.linalg.norm(diff,axis=1)
    idx = np.argpartition(norm, 4)
    tick = 0; ent = 0
    while True:
        if norm[idx[tick]] > 0.01:
            new_coor[ent,0] = coor[idx[tick],0]
            new_coor[ent,1] = coor[idx[tick],1]
            new_coor[ent,2] = coor[idx[tick],2]
            ent += 1     
            #if (coor[idx[tick]] == p_og).all():
                #new_coor[ent,:] = 100
        tick += 1
        if tick == 4:
            break
    return new_coor

#%%
file = open("C168_coordinates", "wb")
# save array to the file
np.save(file, coor)
# close the file
file.close
#%%
      
                
                
        
    




    