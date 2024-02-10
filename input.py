import numpy as np

# Input node matrix, where each row represents the [x,y] coordinates of that node
NODE = np.array([
    [0,0], 
    [0,1], 
    [0,2],
    [2,0],
    [2,1],
    [2,2],
    [4,1],
    [4,2],
    [6,0],
    [6,1],
    [6,2],
    [8,0],
    [8,1],
    [8,2]   
    ])
tnnd = len(NODE) # total number of nodes

# Input element matrix that establishes element-node connectivity
# Each row represents an element: [bottom left, bottom right, top right, top left] 
ELEM = np.array([
    [1,4,5,2], 
    [2,5,6,3], 
    [5,7,8,6],
    [7,10,11,8],
    [10,13,14,11],
    [9,12,13,10]  
    ])
tnel = len(ELEM) # total number of elements

# Young's Modulus of each element
E = np.array([
    10e6,
    10e6,
    10e6,
    10e6,
    20e6,
    20e6
    ])

# Poisson's ratio of each element
Nu = np.array([
    .3,
    .3,
    .3,
    .3,
    .3,
    .3,
    ])

# Thickness of 2D element
t = .1