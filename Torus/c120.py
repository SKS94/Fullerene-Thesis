#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 12:18:25 2020

@author: simon
"""
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

coor = np.zeros([120,3])
coor[:7,:] = np.array([[0.651, -2.917, 0.685], [1.429, -3.125, 1.915], [1.353, -4.575, 2.097], [0, -5.168, 1.792], 
          [0.769, -5.543, -1.548], [1.227, -5.943, -0.214], [0, -5.949, 0.624]])

ax.scatter(coor[:7,0], coor[:7,1], coor[:7,2], s=20, c='b')
count = 7
#%%
reflect_yz = np.identity(3);    reflect_yz[0,0] = -1
num = [0,1,2,4,5]
for ind in num:
    vec = (coor[ind,:]).T 

    new = np.dot(reflect_yz, vec)
    coor[count,:] = new
    ax.scatter(new[0], new[1], new[2], s=20, c='r')
    count += 1
#%%
from numpy import cos, sin
u = np.cos(2*np.pi/5)
v = np.sin(2*np.pi/5)
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
    ax.scatter(xn, yn, zn, s=20, c='y')
    coor[count,0] =  xn;    coor[count,1] =  yn;    coor[count,2] =  zn
    count += 1

#%%
rot_cn_z = array([[np.cos(2*np.pi/5), -np.sin(2*np.pi/5), 0],
            [np.sin(2*np.pi/5), np.cos(2*np.pi/5),  0],
            [0,0,1]])

for i in range(0, count):
    vec = (coor[i,:]).T 
    for j in range(1,5):
        rot_cn_z = array([[np.cos(j*2*np.pi/5), -np.sin(j*2*np.pi/5), 0],
            [np.sin(j*2*np.pi/5), np.cos(j*2*np.pi/5),  0],
            [0,0,1]])
        new = np.dot(rot_cn_z, vec)
        coor[count,:] = new
        ax.scatter(new[0], new[1], new[2], s=20, c='g')
        count+=1
#%%
for i in range(120):
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