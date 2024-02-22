import input # import input module with NODE, ELEM matrices and other parameters
import sympy as sp
import numpy as np

def local_stiffness(elem_num): # for a given element number, will return LK as a matrix
    E = input.E[elem_num-1]
    Nu = input.Nu[elem_num-1]
    element = input.ELEM[elem_num-1] # returns the proper row of the ELEM matrix, which contains the associated nodes
    n1 = element[0] # bottom left node
    n2 = element[1] # bottom right node
    n4 = element[3] # top left node
    a = input.NODE[n2][0] - input.NODE[n1][0] # subtract the x coordinate of n1 from the x coordinate of n2 to get x length of element
    b = input.NODE[n4][1] - input.NODE[n1][1] # subtract the y coordinate of n1 from the y coordinate of n4 to get y length of element
    x = sp.Symbol('x') # define x symbolically
    y = sp.Symbol('y') # define y symbolically
    f1 = (1-(x/a))*(1-(y/b)) # first shape function for 8 DOF rectangular plane stress/plane strain finite element
    f2 = (x/a)*(1-(y/b)) # second shape function for 8 DOF rectangular plane stress/plane strain finite element
    f3 = x*y/a/b # third shape function for 8 DOF rectangular plane stress/plane strain finite element
    f4 = (1-(x/a))*y/b # fourth shape function for 8 DOF rectangular plane stress/plane strain finite element
    shape_function_list = [f1, f2, f3, f4] # create a list so the shape functions can be indexed
    
    def get_shape_functions(i,j): # function used within each of the 3 equations to determine fi and fj and their partial derivatives
        fi = shape_function_list[i-1] # pulls proper shape function for fi
        fj = shape_function_list[j-1] # pulls proper shape function for fj
        dfidx = sp.diff(fi,x) # get partial derivative dfi/dx
        dfjdx = sp.diff(fj,x) # get partial derivative dfj/dx
        dfidy = sp.diff(fi,y) # get partial derivative dfi/dy
        dfjdy = sp.diff(fj,y) # get partial derivative dfj/dy
        return dfidx, dfjdx, dfidy, dfjdy
    
    def double_integral(integrand): # function used to symbolically integrate from 0 to a WRT x, then from 0 to b WRT y
        I1 = sp.integrate(integrand,(x, 0, a)) #integrate with respect to x
        solved_integral = sp.integrate(I1, (y, 0, b)) #integrate with respect to y
        return solved_integral

    kuu = np.zeros((4,4)) # intialize kuu quadrant matrix
    kuv = np.zeros((4,4)) # intialize kuv quadrant matrix
    kvv = np.zeros((4,4)) # intialize kvv quadrant matrix
    multiplier = input.t * E/(1-(Nu**2)) # calculate constant that is multiplied by each stiffness matrix value
    for i in range(1,5): # iterates through each relative index location in the four quadrant matrices of the local stiffness matrix  
        for j in range(1,5):
            dfidx, dfjdx, dfidy, dfjdy = get_shape_functions(i,j) # call function to determine what fi, fj, and their partial derivatives are
            kuu[i-1,j-1] = multiplier * double_integral((dfidx*dfjdx) + (((1-Nu)/2)*dfidy*dfjdy)) # set the given value within the coefficient matrix for kuu
            kuv[i-1,j-1] = multiplier * double_integral((dfidy*dfjdy) + (((1-Nu)/2)*dfidx*dfjdx)) # set the given value within the coefficient matrix for kuv
            kvv[i-1,j-1] = multiplier * double_integral((Nu*dfidx*dfjdy) + (((1-Nu)/2)*dfidy*dfjdx)) # set the given value within the coefficient matrix for kvv
    LK = np.block([[kuu, kuv], [kuv, kvv]]) # Use kuu, kuv, and kvv to form full local stiffness matrix by concatenates quadrant matrices
    return LK
print(local_stiffness(1))