import input as ip
import local_stiffness as lk
import numpy as np

def global_stiffness(): # defines a function to be called to calculate the global stiffness matrix
    GK = np.zeros((2 * ip.tnnd, 2 * ip.tnnd)) # initialize global stiffness matrix
    for k in range(1,ip.tnel+1): # iterates from k = 1:tnel
        print('k: ', k)
        for i in range(1,9): # iterates from i = 1:8 for 2D
            print('i: ', i)
            if i <= 4: # within upper half or Fx portion of the local stiffness matrix
                ii = ip.ELEM[k-1][i-1] # sets the GK row index ii to be the node number of the given LK row 
            else: # within the bottom half or Fy portion of the local stiffness matrix
                i -= 4 # subtracts 4 from i so now Fy1 through Fy4 can be examined
                ii = ip.ELEM[k-1][i-1] + ip.tnnd # adds the total number of elements to account for the shift to the bottom half of GK
            print('ii: ', ii)
            for j in range(1,9): # iterates from j = 1:8 for 2D
                print('j: ', j)
                if j <= 4: # within the left half or u portion of the local stiffness matrix
                    jj = ip.ELEM[k-1][j-1] # sets the GK column index jj to be the node number of the given LK column 
                else: # within the right half or v portion of the local stiffness matrix
                    j -= 4 # subtracts 4 from j so now v1 through v4 can be examined
                    jj = ip.ELEM[k-1][j-1] + ip.tnnd # adds the total number of elements to account for the shift to the right half of GK
                print('jj: ', jj)
                # calls the local_stiffness function within the local_stiffness module, determines LK at element #k, gets LK(i,j), adds it to current value in GK
                GK[ii-1,jj-1] += lk.local_stiffness(k)[i-1,j-1] 
                print(GK)
        
    return GK

GK = global_stiffness()
# print(ip.tnnd)