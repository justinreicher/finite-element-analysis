import input as ip
import local_stiffness as lk
import numpy as np

def global_stiffness(): # defines a function to be called to calculate the global stiffness matrix
    GK = np.zeros((2 * ip.tnnd, 2 * ip.tnnd)) # initialize global stiffness matrix
    for k in range(1,ip.tnel+1): # iterates from k = 1:tnel
        LK = lk.local_stiffness(k)  # precompute local stiffness matrix by calling the local_stiffness function within the local_stiffness module and determining LK at element #k
        element = ip.ELEM[k-1] # predefines the element row
        for i in range(1,9): # iterates from i = 1:8 for 2D
            if i <= 4: # within upper half or Fx portion of the local stiffness matrix
                ii = element[i-1] # sets the GK row index ii to be the node number of the given LK row 
            else: # within the bottom half or Fy portion of the local stiffness matrix
                i -= 4 # subtracts 4 from i so now Fy1 through Fy4 can be examined
                ii = element[i-1] + ip.tnnd # adds the total number of nodes to account for the shift to the bottom half of GK
            for j in range(1,9): # iterates from j = 1:8 for 2D
                if j <= 4: # within the left half or u portion of the local stiffness matrix
                    jj = element[j-1] # sets the GK column index jj to be the node number of the given LK column 
                else: # within the right half or v portion of the local stiffness matrix
                    j -= 4 # subtracts 4 from j so now v1 through v4 can be examined
                    jj = element[j-1] + ip.tnnd # adds the total number of elements to account for the shift to the right half of GK
                GK[ii-1,jj-1] += LK[i-1,j-1] # adds LK(i,j) to current value in GK
    return GK
