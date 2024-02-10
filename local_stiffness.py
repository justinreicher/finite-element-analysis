import input # import input module with NODE, ELEM matrices and other parameters
import sympy as sp
import numpy as np

def local_stiffness(elem_num):
    E = input.E[elem_num-1]
    Nu = input.Nu[elem_num-1]
    element = input.ELEM[elem_num-1] # returns the proper row of the ELEM matrix, which contains the associated nodes
    n1 = element[0] # bottom left node
    n2 = element[1] # bottom right node
    n3 = element[2] # top right node - MIGHT NOT NEED
    n4 = element[3] # top left node
    a = input.NODE[n2][0] - input.NODE[n1][0] # subtract the x coordinate of n1 from the x coordinate of n2 to get x length of element
    b = input.NODE[n4][1] - input.NODE[n1][1] # subtract the y coordinate of n1 from the y coordinate of n4 to get y length of element

    # define x and y symbolically 
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    # define shape functions for 8 DOF rectangular plane stress/plane strain finite element
    f1 = (1-(x/a))*(1-(y/b))
    f2 = (x/a)*(1-(y/b))
    f3 = x*y/a/b
    f4 = (1-(x/a))*y/b
    
    # Calculating kuu, kuv, kvu, kvv
    multiplier = E/(1-Nu**2) # calculate constant that is multiplied by each stiffness matrix value
    def get_shape_functions(i,j): # function used within each of the 3 equations to determine fi and fj 
        if i == 1:
                fi = f1
        elif i == 2:
            fi = f2
        elif i == 3:
            fi = f3
        else:
            fi = f4
        if j == 1:
            fj = f1
        elif j == 2:
            fj = f2
        elif j == 3:
            fj = f3
        else:
            fj = f4
        
        # get partial derivatives
        dfidx = sp.diff(fi,x)
        dfjdx = sp.diff(fj,x)
        dfidy = sp.diff(fi,y)
        dfjdy = sp.diff(fj,x)
        
        return fi, fj, dfidx, dfjdx, dfidy, dfjdy
    
    kuu = np.zeros((4,4))
    kuv = kuu
    kvv = kuu
    for i in range(1,5):
        for j in range(1,5):
            fi, fj, dfidx, dfjdx, dfidy, dfjdy = get_shape_functions(i,j) # call function to determine what fi, fj, and their partial derivatives are
            
            integranduu = ((dfidx*dfjdx) + (((1-Nu)/2)*dfidy*dfjdy)) * input.t
            I1uu = sp.integrate(integranduu,(x, 0, a)) #integrate with respect to x
            I2uu = sp.integrate(I1uu, (y, 0, b)) #integrate with respect to y
            kuu[i-1,j-1] = multiplier * I2uu # set the given value within the coefficient matrix

            integranduv = ((dfidy*dfjdy) + (((1-Nu)/2)*dfidx*dfjdx)) * input.t
            I1uv = sp.integrate(integranduv,(x, 0, a)) #integrate with respect to x
            I2uv = sp.integrate(I1uv, (y, 0, b)) #integrate with respect to y
            kuv[i-1,j-1] = multiplier * I2uv # set the given value within the coefficient matrix

            integrandvv = ((Nu*dfidx*dfjdy) + (((1-Nu)/2)*dfidy*dfjdx)) * input.t
            I1vv = sp.integrate(integrandvv,(x, 0, a)) #integrate with respect to x
            I2vv = sp.integrate(I1vv, (y, 0, b)) #integrate with respect to y
            kvv[i-1,j-1] = multiplier * I2vv # set the given value within the coefficient matrix

    # Use kuu, kuv, and kvv to form full local stiffness matrix
    LK = np.vstack((np.hstack((kuu, kuv)),np.hstack((kuv, kvv)))) # concatenates quadrant matrices
    print(kuu)
    print(kuv)
    print(kvv)
    return LK
            

LK = local_stiffness(1)
print(LK)
print(LK[0,7])