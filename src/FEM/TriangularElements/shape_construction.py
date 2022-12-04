#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:30:59 2020

@author: simon
"""

def FE_construction(shape):
    import numpy as np
    if shape == 'triangle_equil':
        Element = np.array([[-1/2, np.sqrt(3)/2], [1, 0], [-1/2, -np.sqrt(3)/2]])
        determinant = np.sqrt(27)/2
    elif shape == 'triangle_isosceles':
        Element = np.array([[-1/2, np.sqrt(3)/2], [1/2 *np.sqrt(3)/np.sin(72/2*np.pi/180)*np.sin(54*np.pi/180) - 1/2, 0], [-1/2, -np.sqrt(3)/2]])
        determinant = 2*(np.sqrt(3)/2 * 1/2 *np.sqrt(3)/np.sin(72/2*np.pi/180)*np.sin(54*np.pi/180)) #2*area
        #Element = np.roll(Element,1,axis=0)
    elif shape == 'triangle_isoceles_heptagon':
        Element = np.array([[-1/2, np.sqrt(3)/2], [1.7983202809335845-1/2, 0], [-1/2, -np.sqrt(3)/2]])
        determinant = 2 * np.sqrt(3)/2 * 1.7983202809335845  #2*area
    
    '''If full structure consits of several different elements(quite likely), then several elements
    should be returned #change.'''
    return determinant, Element


'''

import numpy as np
from matplotlib import pyplot as plt
a = 72*np.pi/180
det, Element = FE_construction('triangle_isosceles')
det, Ei = FE_construction('triangle_isosceles')
det, E = FE_construction('triangle_isosceles')
det, Ee = FE_construction('triangle_equil')
#R = np.array([[np.cos(a), -np.sin(a)],[np.sin(a), np.cos(a)]])

plt.scatter(E[:,0], E[:,1])
#plt.scatter([x_l - 0.5959908542188245, x_l -  0.5959908542188245], [0.43301270189221924, - 0.43301270189221924])
#plt.scatter((-E[0,0] + E[1,0])/2, E[0,1]/2)
#plt.scatter((E[0,0] + E[1,0])/2, -E[0,1]/2)

plt.scatter((-E[0,0] + E[1,0])/2-0.5, -E[0,1]/2)
plt.scatter((-E[0,0] + E[1,0])/2-0.5, E[0,1]/2)

plt.plot(E[:,0], E[:,1])

def tester(A, *argv):
    print(A)
    print(len(argv))
    #for arg in argv:
        #print(arg)
    print(argv[0])
    print(argv[0][0,0])
        
#tester(1,3,4,5,6,6)
#tester(1)

tester(A,B)
'''
'''
def rot(coor, theta): #rotation, used for reflection as well.
    rot = np.asarray([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    return np.dot(rot, coor)


plt.scatter(E[:,0], E[:,1])
for i in range(3):
    E[i,:] = rot(E[i,:], a)
plt.scatter(E[:,0], E[:,1])


for i in range(3):
    E[i,:] = rot(E[i,:], a)
plt.scatter(E[:,0], E[:,1])


for i in range(3):
    E[i,:] = rot(E[i,:], a)
plt.scatter(E[:,0], E[:,1]) 

for i in range(3):
    E[i,:] = rot(E[i,:], a)
plt.scatter(E[:,0], E[:,1])
'''

#%%
'''
coor = Ei[2,:]
trans_mat = np.array([[Ei[1,0] - Ei[0,0], Ei[2,0] - Ei[0,0]], 
                          [Ei[1,1] - Ei[0,1], Ei[2,1] - Ei[0,1]]])
    
inv_Jacobian = np.linalg.inv(trans_mat)
quad_ref_sys = np.dot(inv_Jacobian,   coor - Element[0,:]) #Quadrature points in canonical system

coordinates_7 = np.array([[ 0.,  0.],        [-0.41042619,  0.        ], [ 0.2052131 ,  0.35543951], 
       [ 0.2052131 , -0.35543951], [ 0.69614048,  0.], [-0.34807024, -0.60287534], [-0.34807024,  0.60287534]])
weights_7 = np.array([0.29228357, 0.17198505, 0.17198505, 0.17198505, 0.16359979, 0.16359979, 0.16359979])

plt.scatter(coordinates_7[:,0], coordinates_7[:,1])

plt.scatter(coordinates_7[:,0] + (1/2 *np.sqrt(3)/np.sin(72/2*np.pi/180)*np.sin(54*np.pi/180)-1), coordinates_7[:,1])
#%%
Element = FE_construction('triangle_equil')
trans_mat = np.array([[Element[1,0] - Element[0,0], Element[2,0] - Element[0,0]], 
                          [Element[1,1] - Element[0,1], Element[2,1] - Element[0,1]]])

'''


'''
import numpy as np
det, Element = FE_construction('triangle_isosceles')
det1, Element1 = FE_construction('triangle_equil')
trans_mat = np.array([[Element[1,0] - Element[0,0], Element[2,0] - Element[0,0]], 
                          [Element[1,1] - Element[0,1], Element[2,1] - Element[0,1]]])

    
print(Element1[0,:] + np.dot(trans_mat, np.array([[1/3],[1/3]])).T)

from matplotlib import pyplot as plt
plt.scatter(Element[:,0], Element[:,1],c='b')


def A(coor):
    det, Element = FE_construction('triangle_equil')
    trans_mat = np.array([[Element[1,0] - Element[0,0], Element[2,0] - Element[0,0]], 
                      [Element[1,1] - Element[0,1], Element[2,1] - Element[0,1]]])

    inv_Jacobian = np.linalg.inv(trans_mat)
    quad_ref_sys = np.dot(inv_Jacobian,   coor - Element[0,:]) #Quadrature points in canonical system
    return quad_ref_sys

for i in range(7):
    
    a = A(coordinates_7[i,:])
    
    b = Element1[0,:] + np.dot(trans_mat, a)
    plt.scatter(b[0], b[1], c='b')
#%%
    
a_1 = np.zeros([3,3])    
a_1[:2,:] = Element1.T
a_1[2,:] = 1

a_2 = np.zeros([3,3])    
a_2[:2,:] = Element.T
a_2[2,:] = 1    

T = np.dot(a_2, np.linalg.inv(a_1))

coor = np.array([-0.5,0,1])
print(np.dot(T,coor))

plt.figure(3)
for i in range(7):
    coor = np.array([coordinates_7[i,0], coordinates_7[i,1],1])
    res = np.dot(T,coor)
    
    trans_mat = np.array([[Element[1,0] - Element[0,0], Element[2,0] - Element[0,0]], 
                          [Element[1,1] - Element[0,1], Element[2,1] - Element[0,1]]])
    inv_J = np.linalg.inv(trans_mat)
    res = np.dot(inv_J,   res[:2] - Element[0,:])
    plt.scatter(res[0], res[1], c='g')
    
    
    trans_mat = np.array([[Element1[1,0] - Element1[0,0], Element1[2,0] - Element1[0,0]], 
                          [Element1[1,1] - Element1[0,1], Element1[2,1] - Element1[0,1]]])
    inv_J = np.linalg.inv(trans_mat)
    res = np.dot(inv_J,   coordinates_7[i,:] - Element1[0,:])
    plt.scatter(res[0], res[1], c='b')
#%%
for i in range(3):
    #coor = np.ones(3)
    coor = Element[i,:]
    
#    res = np.dot(T,coor)
    
    trans_mat = np.array([[Element[1,0] - Element[0,0], Element[2,0] - Element[0,0]], 
                          [Element[1,1] - Element[0,1], Element[2,1] - Element[0,1]]])
    inv_J = np.linalg.inv(trans_mat)
    res = np.dot(inv_J,   coor - Element[0,:])
    print('At point ' + str(res))
    print(Basis.BasisFunctionValue(5, res[:2]))
    
    
'''
    
    