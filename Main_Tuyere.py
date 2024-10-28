import numpy as np
from Mesh_Class import Mesh
from Solver_Class import Solver
import matplotlib.pyplot as plt


# Constantes
gamma = 1.4
x1 = 0.0  
x2 = 10.0
n = 201
t = 250.0  
CFL = 1.1
method = "Macormack"

#Cas tube: Aire constante
def area_tuyere(x): return 1.398 + 0.347*np.tanh(0.8 * x - 4)

Aire = []
for i in np.linspace(x1,x2,n):
    Aire.append(area_tuyere(i))


# Conditions amont (entrée)
rho_0 = 1.25 # Densité amont
p_0 = 1.0    # Pression amont
gamma = 1.4
a_0 = np.sqrt(gamma * p_0 / rho_0)  # Vitesse du son amont
M_inlet = 1.25
u_0 = M_inlet * a_0  # Vitesse à l'entrée

#Conditions initiale
x = np.linspace(x1,x2,n)
init_p = np.zeros(n)
init_rho = np.zeros(n)
init_u = np.zeros(n)
init_e = np.zeros(n)

init_rho = np.full(n, rho_0)
init_p = np.full(n, p_0)
init_u = np.full(n, u_0)
init_e = init_p / (gamma - 1) + 0.5 * init_rho * init_u**2







Maillage = Mesh(x1, x2, n, Aire, t)


Solveur = Solver(Maillage,method,CFL)
Solveur.Conditions_Initiales(init_p,init_rho,init_u,init_e)
Solveur.Solve_Macormack()






    
    

