"""module for computing the global stiffness matrix"""

import input as ip
import local_stiffness as lk
import numpy as np
import pandas as pd

def global_stiffness(): # defines a function to be called to calculate the global stiffness matrix
    """
    calculate the global stiffness matrix

    returns:
        numpy.ndarray: the global stiffness matrix
    """
    # initialize global stiffness matrix
    GK = np.zeros((2 * ip.tnnd, 2 * ip.tnnd))

    # iterate through each element and store its contribution to the global stiffness matrix
    for k in range(1,ip.tnel+1):
        # precompute the local stiffness matrix
        LK = lk.local_stiffness(k)  

        # determine element node connectivity
        element = ip.ELEM[k-1]

        # iterate through all 8 DOF of the element
        for i in range(1,9):
            # determine the row index (ii) in the global stiffness matrix
            if i <= 4: # within upper half or Fx portion of the local stiffness matrix
                ii = element[i-1] # sets the GK row index ii to be the node number of the given LK row 
            else: # within the bottom half or Fy portion of the local stiffness matrix
                i -= 4 # subtracts 4 from i so now Fy1 through Fy4 can be examined
                ii = element[i-1] + ip.tnnd # adds the total number of nodes to account for the shift to the bottom half of GK
            
            # iterate through all 8 DOF of the element
            for j in range(1,9):
                # determine the column index (jj) in the global stiffness matrix
                if j <= 4: # within the left half or u portion of the local stiffness matrix
                    jj = element[j-1] # sets the GK column index jj to be the node number of the given LK column 
                else: # within the right half or v portion of the local stiffness matrix
                    j -= 4 # subtracts 4 from j so now v1 through v4 can be examined
                    jj = element[j-1] + ip.tnnd # adds the total number of elements to account for the shift to the right half of GK
                
                # add the contribution of the local stiffness matrix to the current value in the global stiffness matrix at the computed indeces
                GK[ii-1,jj-1] += LK[i-1,j-1] # adds LK(i,j) to current value in GK
    
    return GK

if __name__ == "__main__":
    """
    will only run when the script is executed directly, not when it's imported as a module in another script
    """
    # calculate the global stiffness matrix
    GK = global_stiffness()
    
    # print the global stiffness matrix
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    print("Global stiffness matrix (GK):\n", pd.DataFrame(GK))