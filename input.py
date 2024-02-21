import numpy as np

# Input node matrix, where each row represents the [x,y] coordinates of that node
NODE = np.array([
    [0,0], 
    [0,2], 
    [0,4],
    [3,0],
    [3,2],
    [3,4],
    [6,2],
    [6,4],
    [9,2],
    [9,4],
    [12,0],
    [12,2],
    [12,4],
    [15,0],
    [15,2],
    [15,4] 
    ])
tnnd = len(NODE) # total number of nodes

# Input element matrix that establishes element-node connectivity
# Each row represents an element: [bottom left, bottom right, top right, top left] 
ELEM = np.array([
    [1,4,5,2], 
    [2,5,6,3], 
    [5,7,8,6],
    [7,9,10,8],
    [9,12,13,10],
    [12,15,16,13],
    [11,14,15,12]  
    ])
tnel = len(ELEM) # total number of elements

# Young's Modulus of each element
E = np.array([
    10e6,
    10e6,
    10e6,
    10e6,
    10e6,
    10e6,
    10e6
    ])

# Poisson's ratio of each element
Nu = np.array([
    .3,
    .3,
    .3,
    .3,
    .3,
    .3,
    .3
    ])

# Thickness of 2D element
t = .1