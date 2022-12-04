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