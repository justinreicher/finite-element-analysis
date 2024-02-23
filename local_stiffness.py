"""module for computing local stiffness matrices"""

import input # import input module with NODE, ELEM matrices and other parameters
import sympy as sp
import numpy as np

def local_stiffness(elem_num):
    """
    calculate the local stiffness matrix for a given element number.

    args:
        elem_num (int): the element number of interest
    
    returns:
        numpy.ndarray: the local stiffness matrix for the element of interest
    """
    # material properites and element node connectivity
    E = input.E[elem_num-1]
    Nu = input.Nu[elem_num-1]
    element = input.ELEM[elem_num-1] # returns the proper row of the ELEM matrix, which contains the associated nodes
    n1, n2, _, n4 = element # saves the bottom left, bottom right, and top left node numbers

    # calculate element dimensions
    a = input.NODE[n2][0] - input.NODE[n1][0] # subtract the x coordinate of n1 from the x coordinate of n2 to get x length of element
    b = input.NODE[n4][1] - input.NODE[n1][1] # subtract the y coordinate of n1 from the y coordinate of n4 to get y length of element

    # initialize symbolic variables x and y
    x = sp.Symbol('x')
    y = sp.Symbol('y')

    # define shape functions for 8 DOF rectangular plane stress/strain finite element
    f1 = (1 - (x / a)) * (1 - (y / b))
    f2 = (x / a) * (1 - (y / b))
    f3 = x * y / (a * b)
    f4 = (1 - (x / a)) * (y / b)
    shape_function_list = [f1, f2, f3, f4]
    
    def get_shape_functions(i,j):
        """return the shape function's partial derivatives"""
        fi = shape_function_list[i-1]
        fj = shape_function_list[j-1]
        dfidx = sp.diff(fi,x) # dfi/dx
        dfjdx = sp.diff(fj,x) # dfj/dx
        dfidy = sp.diff(fi,y) # dfi/dy
        dfjdy = sp.diff(fj,y) # dfj/dy
        return dfidx, dfjdx, dfidy, dfjdy
    
    def double_integral(integrand):
        """symbolic double integration over the area of the 2D element"""
        I1 = sp.integrate(integrand,(x, 0, a)) #integrate with respect to x
        solved_integral = sp.integrate(I1, (y, 0, b)) #integrate with respect to y
        return solved_integral

    # initialize quadrant kuu, kuv, and kvv matrices
    kuu = np.zeros((4,4))
    kuv = np.zeros((4,4))
    kvv = np.zeros((4,4))

    # calculate constant that is pulled out of the area integrals
    multiplier = input.t * E/(1-(Nu**2))

    # compute the local stiffness matrix by iterating through each relative index location in the four quadrant matrices of the local stiffness matrix  
    for i in range(1,5):
        for j in range(1,5):
            dfidx, dfjdx, dfidy, dfjdy = get_shape_functions(i,j) # determine shape function partial derivatives
            kuu[i-1,j-1] = multiplier * double_integral((dfidx * dfjdx) + (((1 - Nu) / 2) * dfidy * dfjdy))
            kvv[i-1,j-1] = multiplier * double_integral((dfidy * dfjdy) + (((1 - Nu) / 2) * dfidx * dfjdx))
            kuv[i-1,j-1] = multiplier * double_integral((Nu * dfidx * dfjdy) + (((1 - Nu) / 2) * dfidy * dfjdx))

    # form full local stiffness matrix by concatenating quadrant matrices
    LK = np.block([[kuu, kuv], [kuv, kvv]]) 
    return LK


if __name__ == "__main__":
    """
    will only run when the script is executed directly, not when it's imported as a module in global_stiffness.py
    """
    # test the local_stiffness function by calculating the local stiffness matrix for element #1
    LK = local_stiffness(1)

    # print the local stiffness matrix for element #1
    print("Local Stiffness Matrix (LK) for element #1:\n", LK)