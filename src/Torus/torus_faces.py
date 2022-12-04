#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 11:29:33 2020

@author: simon
"""
import numpy as np
file = open("C168_coordinates", "rb")
coor = np.load(file)
file.close()
from matplotlib import pyplot as plt; from mpl_toolkits.mplot3d import Axes3D
from numpy import array
plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#%%
def neighbour(coor, point):
    new_coor = np.zeros([3,3])
    diff = (coor - point)
    norm = np.linalg.norm(diff,axis=1)
    idx  = np.argpartition(norm, 4)
    tick = 0; ent = 0
    indices = []
    while True:
        if norm[idx[tick]] > 0.01:
            new_coor[ent,0] = coor[idx[tick],0]
            new_coor[ent,1] = coor[idx[tick],1]
            new_coor[ent,2] = coor[idx[tick],2]
            indices.append(idx[tick])
            ent += 1     
            #if (coor[idx[tick]] == p_og).all():
                #new_coor[ent,:] = 100
        tick += 1
        if tick == 4:
            break
    return new_coor, indices
#%%
'''
for i in range(1):
    point = coor[i,:]
    coor_new = neighbour(coor, point)
    prev_point = point
    for j in range(3):
        while True:
            c = coor_new[j]
            c_new1 = neighbour(coor, c)
            c_new = np.zeros([2,3]); tick = 0
            for temp in c_new1:
                if ((temp == prev_point).all()) == False:
                    c_new[tick,:] = temp
                    tick+=1
            break
'''
#%%
'''
p = point
cp = c - p
a = c_new[0,:]; b = c_new[1,:]

ax.plot([p[0], c[0]],[p[1], c[1]],[p[2], c[2]], c='b')

ax.plot([a[0], c[0]],[a[1], c[1]],[a[2], c[2]], c='r')
ax.plot([b[0], c[0]],[b[1], c[1]],[b[2], c[2]], c='g')

a_red = a - c;  b_red = b - c;  

angle_a = np.arccos(np.dot(cp,a_red)/(np.linalg.norm(cp)*np.linalg.norm(a_red)))
angle_b = np.arccos(np.dot(cp,b_red)/(np.linalg.norm(cp)*np.linalg.norm(b_red)))

print()'''
#%%
'''
import math
s = np.linalg.norm(np.cross(cp,aa))
c = np.dot(cp,aa)
angle = math.atan2(s, c)
print(angle)

s = np.linalg.norm(np.cross(cp,bb))
c = np.dot(cp,bb)
angle = math.atan2(s, c)
print(angle)
'''
#%% Mean-point method
'''
mean_p = (a+b)/2
n = np.cross(cp, mean_p)
res = np.dot(n, b_red)
print(res)
'''
faces = []; faces_plot = []
for i in range(168):
    point_og = coor[i,:]
    if i == 2:
            print('a')
    coor_new, indices1 = neighbour(coor, point_og)
    for j in range(3):
        
        c = coor_new[j,:]
        prev_point = np.copy(point_og)
        cul = c + prev_point
        big_p = np.zeros([7,3]); k = 0
        
        face_ind = [i, indices1[j]]
        while True:
            pos_indices = []
            c_neighbours, indices = neighbour(coor, c)
            c_neighbours_mod = np.zeros([2,3]); tick = 0
            for h in range(3):
                if ((c_neighbours[h] == prev_point).all()) == False:
                    c_neighbours_mod[tick,:] = c_neighbours[h]
                    pos_indices.append(indices[h])
                    
                    tick+=1
            
            a = c_neighbours_mod[0,:];  b = c_neighbours_mod[1,:]
            mean_point = (a+b)/2
            
            n = np.cross(c - prev_point, mean_point - prev_point)
            res_a = np.dot(n, a - prev_point)
            res_b = np.dot(n, b - prev_point)
            
            big_p[k, :] = prev_point; k+=1
            #if k == 7:
             #break
            prev_point = np.copy(c)
            if res_a < 0:
                c = np.copy(a); face_ind.append(pos_indices[0])
            else:
                c = np.copy(b); face_ind.append(pos_indices[1])
            
            if face_ind[0] == face_ind[-1]:
                #face_ind.pop()
                break
            
        faces.append(face_ind[:-1]); faces_plot.append(face_ind)
        ax.plot(coor[face_ind,0], coor[face_ind,1], coor[face_ind,2])        
#%%

ind = 0
while True:
    try:
        max_iter = len(faces) - ind;    i = 1
        #for i in range(1, len(faces) - ind):
        while True:
            if ind == 81 and i == 30:
                print('here')
            if set(faces[ind]).issubset(faces[-i]) and faces[ind] != faces[-i]:
                faces.pop(-i)
                #max_iter -= 1
                max_iter = len(faces) - ind
            else:
                i+=1
            if i > max_iter:
                break
    except IndexError:
        break
    ind += 1
#%%
'''
for i in range(len(faces)):
    if len(faces[i]) == 7:
        ax.plot(coor[faces[i],0], coor[faces[i],1], coor[faces[i],2])
        print(i)
        print(faces[i])
'''

def neighbour_alt(coor, point, num):
    new_coor = np.zeros([num,3])
    diff = (coor - point)
    norm = np.linalg.norm(diff,axis=1)
    idx  = np.argpartition(norm, num)
    tick = 0; ent = 0
    indices = []
    while True:
        if norm[idx[tick]] > 0.01:
            new_coor[ent,0] = coor[idx[tick],0]
            new_coor[ent,1] = coor[idx[tick],1]
            new_coor[ent,2] = coor[idx[tick],2]
            indices.append(idx[tick])
            ent += 1     
            #if (coor[idx[tick]] == p_og).all():
                #new_coor[ent,:] = 100
        if tick == num:
            break
        tick += 1
    return new_coor, indices

#%%
dual_coor = np.zeros([len(faces), 3]); i = 0
for face in faces:
    face_coor = coor[face]
    dual_coor[i,:] = np.mean(face_coor, axis=0)
    i+=1
ax.scatter(dual_coor[:,0], dual_coor[:,1], dual_coor[:,2], c='r')
#%%
dual_faces = []; dual_faces_plot = []
for i in range(dual_coor.shape[0]):
    #neighbours, indices = neighbour_alt(dual_coor, dual_coor[i,:], len(faces[i]))
    indices = faces[i]
    for j in range(len(faces[i])):
        pairs = np.roll(indices,j)[0:3]
        face_ind = [i];  k = 0
        while True:
            if k!=i:
                if pairs[0] in faces[k] and pairs[1] in faces[k]:
                    face_ind.append(k)
                elif pairs[1] in faces[k] and pairs[2] in faces[k]:
                    face_ind.append(k)    
            if len(face_ind)==3:
                break    
            k+=1
        dual_faces.append(face_ind);   dual_faces_plot.append(face_ind)
#%%
ind = 0
while True:
    try:
        max_iter = len(dual_faces) - ind;    i = 1
        #for i in range(1, len(faces) - ind):
        while True:
            if set(dual_faces[ind]).issubset(dual_faces[-i]) and dual_faces[ind] != dual_faces[-i]:
                dual_faces.pop(-i)
                #max_iter -= 1
                max_iter = len(dual_faces) - ind
            else:
                i+=1
            if i > max_iter:
                break
    except IndexError:
        break
    ind += 1
#%%
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
fig = plt.figure(figsize=plt.figaspect(0.5)*1.5)
ax = fig.add_subplot(111, projection='3d')
cmap = plt.get_cmap("tab10")
for face in dual_faces_plot:
    #ax.plot(dual_coor[face][:,0], dual_coor[face][:,1], dual_coor[face][:,2], c='r')
    x = dual_coor[face][:,0]
    y = dual_coor[face][:,1]
    z = dual_coor[face][:,2]
    verts = [list(zip(x,y,z))]
    ax.add_collection3d(Poly3DCollection(verts, alpha = 1, edgecolors='lightgrey', facecolor=cmap(2)))
    
    print(dual_coor[face])
ax.set_xlim([-8,8])
ax.set_ylim([-8,8])
ax.set_zlim([-2,2])
ax.axis('off')
#ax.set_aspect('equal')
ax.view_init(elev=50., azim=60)

#%%
for i in range(5):
    c = dual_coor[dual_faces[i]]
    print(np.linalg.norm(c[0,:] - c[1,:]))
    print(np.linalg.norm(c[1,:] - c[2,:]))

#%%
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
fig = plt.figure(figsize=plt.figaspect(0.5)*1.5)
ax = fig.add_subplot(111, projection='3d')
a = 0.9
cmap = plt.get_cmap("tab10")
pent = 0
hep = 0
for i in range(len(faces)):
    x = coor[faces[i],0]
    y = coor[faces[i],1]
    z = coor[faces[i],2]
    verts = [list(zip(x,y,z))]
    if len(faces[i]) == 5:
        ax.add_collection3d(Poly3DCollection(verts, alpha = a, edgecolors='lightgrey', facecolor=cmap(0)))
        pent+=1
    elif len(faces[i]) == 6:
        ax.add_collection3d(Poly3DCollection(verts, alpha = a, edgecolors='lightgrey', facecolor=cmap(1)))
    elif len(faces[i]) == 7:
        ax.add_collection3d(Poly3DCollection(verts, alpha = a, edgecolors='lightgrey', facecolor=cmap(3)))
        hep+=1
    
    
    #ax.scatter(np.mean(x), np.mean(y), np.mean(z), c='orange', s= 40)
ax.set_xlim([-8,8])
ax.set_ylim([-8,8])
ax.set_zlim([-2,2])


plt.axis('off')