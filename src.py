"""main module that serves as the entry point for conducting FEA analysis"""

import global_stiffness as gs

# compute global stiffness matrix
GK = gs.global_stiffness()
