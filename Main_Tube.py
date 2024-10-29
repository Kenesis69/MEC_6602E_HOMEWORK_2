import numpy as np
from Mesh_Class import Mesh
from Solver_Class import Solver
import matplotlib.pyplot as plt


# Constantes
gamma = 1.4
x1 = 0.0  
x2 = 1000.0  
n = 1001
t = 250.0  
CFL = 0.9
method = "Macormack"

#Cas tube: Aire constante
Aire = np.linspace(1,1,n)
#Conditions initiale
x = np.linspace(x1,x2,n)
init_p = np.zeros(n)
init_rho = np.zeros(n)
init_u = np.zeros(n)
init_e = np.zeros(n)

for i in range(0,n):
    if x[i] < 500:
        init_p[i] = 1
        init_rho[i] = 1
        init_u[i] = 0
        init_e[i] = init_p[i]/(gamma -1) + init_rho[i]*init_u[i]**2 * 0.5
    elif x[i] >= 500: 
        init_p[i] = 4
        init_rho[i] = 4
        init_u[i] = 0
        init_e[i] = init_p[i]/(gamma -1) + init_rho[i]*init_u[i]**2 * 0.5







Maillage = Mesh(x1, x2, n, Aire, t)


Solveur = Solver(Maillage,method,CFL)
Solveur.Conditions_Initiales(init_p,init_rho,init_u,init_e)
Solveur.Solve_Macormack()






    
    

