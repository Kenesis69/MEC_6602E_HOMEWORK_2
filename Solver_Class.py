import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

import Mesh_Class as Mesh


class Solver():

    def __init__(self, mesh, method, CFL):
        self.mesh = mesh
        self.method = method
        self.CFL = CFL 
    

    def Conditions_Initiales(self,pression,rho,u,e,gamma = 1.4):
        self.p = pression
        self.rho = rho
        self.u = u
        self.e = e 
        self.gamma = gamma
        print("================ Here are the initials conditions ============================")
        print("==============================================================================")
        print("============================== Pressure ======================================")
        print (pression)
        print("============================== Rho ===========================================")
        print (rho)
        print("============================== Speed =========================================")
        print (u)
        print("============================== Energy ========================================")
        print (e)
    


    def Solve_Macormack(self):
        # We have to initialise Q,F and S
        Q = np.zeros((3,self.mesh.n))
        F = np.zeros((3,self.mesh.n))
        S = np.zeros((3,self.mesh.n))
        Q_predictor = np.zeros((3,self.mesh.n+2)) #Already cell Ghosted
        F_predictor = np.zeros((3,self.mesh.n+2)) #Already cell Ghosted
        S_predictor = np.zeros((3,self.mesh.n+2)) #Already cell Ghosted
        Q_corrector = np.zeros((3,self.mesh.n+2)) #Already cell Ghosted
        F_corrector = np.zeros((3,self.mesh.n+2)) #Already cell Ghosted
        S_corrector = np.zeros((3,self.mesh.n+2)) #Already cell Ghosted

        #Formulas for Q, F, and S
        Q[0] = self.rho*self.mesh.area
        Q[1] = self.rho*self.u*self.mesh.area
        Q[2] = self.e*self.mesh.area
        F[0] = Q[1]
        F[1] = (self.rho * self.u**2 + self.p)*self.mesh.area
        F[2] = self.u * (self.e + self.p) * self.mesh.area
        S[0] = 0
        S[1] = self.p*self.mesh.dAdx
        S[2] = 0

        #Cell Ghosting
        Left_Col_Q = Q[:,[0]]
        Right_Col_Q = Q[:,[-1]]
        Left_Col_F= F[:,[0]]
        Right_Col_F = F[:,[-1]]
        Left_Col_S = S[:,[0]]
        Right_Col_S = S[:,[-1]]
        Q = np.hstack((Left_Col_Q,Q,Right_Col_Q))
        F = np.hstack((Left_Col_F,F,Right_Col_F))
        S = np.hstack((Left_Col_S,S,Right_Col_S))
        init_Q = Q


        print("================ Here are the initial Vectors ============================")
        print("==============================================================================")
        print("============================== Q ======================================")
        print (Q)
        print("============================== F ===========================================")
        print (F)
        print("============================== S =========================================")
        print (S)
        print("============================== SIZE =========================================")
        print("Cell Ghosted",Q.shape[1])
        print("Non cell ghosted",Q.shape[1]-2)

        # What is the first time space acoording to the CFL ? 
        Delta_t = self.CFL/(abs(self.u) + (self.gamma*self.p/self.rho)**0.5)*self.mesh.delta_x
        delta_t = min(Delta_t)
        live_time = 0
        iteration = 0
        iterations = []
        L2s = []   
        # Now the begin the Solver
        
        # Activer le mode interactif
            # Activer le mode interactif
        plt.ion()
        fig, (ax1, ax2, ax3,ax4) = plt.subplots(4, 1, figsize=(8, 8))
        fig.suptitle("distribution of density, pressure ,mach number and L2, along x")
        # Initialiser les lignes pour la densité et la pression
        line1, = ax1.plot(self.mesh.x, self.rho, label='Densité ')
        ax1.set_xlabel('x')
        ax1.set_ylabel('Densité')
        ax1.legend()
        ax1.grid(True)

        line2, = ax2.plot(self.mesh.x, self.p, label='Pression ')
        ax2.set_xlabel('x')
        ax2.set_ylabel('Pression ')
        ax2.legend()
        ax2.grid(True)


        line3, = ax3.plot(self.mesh.x, self.u/(self.gamma*self.p/self.rho)**0.5, label='Mac Number ')
        ax3.set_xlabel('x')
        ax3.set_ylabel('Mach number')
        ax3.legend()
        ax3.grid(True)

        line4, = ax4.plot([], [], label='L2 Error')
        ax4.set_xlabel('log Iterations')
        ax4.set_ylabel('log L2 Error')
        ax4.legend()
        ax4.grid(True)



        plt.tight_layout()

        while live_time < self.mesh.t:
            print("=============================== ITERATION, ", iteration," ===========================")
            print("=============================== Delta_t, ", delta_t," ===============================")
            print("=============================== Live Time, ", live_time," ===============================")
            #predictor step
            #Q_predictor[:,1:-1] = Q[:,1:-1] - delta_t*(F[:,2:] - F[:,1:-1])/self.mesh.delta_x + delta_t * S[:,1:-1]
            Q_predictor[:,1:-1] = Q[:,1:-1] - delta_t*(F[:,1:-1] - F[:,0:-2])/self.mesh.delta_x + delta_t * S[:,1:-1]
            # Calcul de la viscosité artificielle
            
            #properties on predictor stage
            rho_predictor = Q_predictor[0,1:-1]/self.mesh.area
            u_predictor = Q_predictor[1,1:-1]/(rho_predictor*self.mesh.area)
            e_predictor = Q_predictor[2,1:-1]/self.mesh.area
            p_predictor = (self.gamma - 1)*(e_predictor - 0.5*(rho_predictor*u_predictor**2))
            
            

            # F_predictor S_predictor
            F_predictor[0,1:-1] = rho_predictor*u_predictor*self.mesh.area
            F_predictor[1,1:-1] = (rho_predictor*u_predictor**2 +  p_predictor)*self.mesh.area
            F_predictor[2,1:-1] = u_predictor * ( e_predictor+p_predictor) * self.mesh.area
            S_predictor[1,1:-1] = p_predictor * self.mesh.dAdx
            F_predictor[:,0] = F_predictor[:,1]
            F_predictor[:,-1] = F_predictor[:,-2]
            S_predictor[:,0] = S_predictor[:,1]
            S_predictor[:,-1] = S_predictor[:,-2]


            # Corrector step
            #Q_corrector[:, 1:-1] = Q[:, 1:-1] - delta_t * (F_predictor[:, 1:-1] - F_predictor[:, 0:-2]) / self.mesh.delta_x + delta_t * S_predictor[:, 1:-1]
            Q_corrector[:, 1:-1] = Q[:, 1:-1] - delta_t * (F_predictor[:, 2:] - F_predictor[:, 1:-1]) / self.mesh.delta_x + delta_t * S_predictor[:, 1:-1]
            Q_corrector = (Q_predictor + Q_corrector)*0.5


            # Re-implementation conditions LIMITES au cellghost
            #Q_corrector[:,0] = Q[:,1]  # TUBE 
            Q_corrector[:,0] = init_Q[:,0]  # super sonic in 
            Q_corrector[:,-1] = Q[:,-2] # TUBE ET SUPERSONIC OUT

            #properties on predictor stage
            rho_corrector = Q_corrector[0,1:-1]/self.mesh.area
            u_corrector = Q_corrector[1,1:-1]/(rho_predictor*self.mesh.area)
            e_corrector = Q_corrector[2,1:-1]/self.mesh.area
            p_corrector = (self.gamma - 1)*(e_corrector - 0.5*(rho_corrector*u_corrector**2))
            
            
            
            


            # F_corrector S_corrector
            F_corrector[0,1:-1] = rho_corrector*u_corrector*self.mesh.area
            F_corrector[1,1:-1] = (rho_corrector*u_corrector**2 +  p_corrector)*self.mesh.area
            F_corrector[2,1:-1] = u_corrector * ( e_corrector+p_corrector) * self.mesh.area
            S_corrector[0,1:-1] = 0
            S_corrector[1,1:-1] = p_corrector * self.mesh.dAdx
            S_corrector[2,1:-1] = 0
            F_corrector[:,0] = F_corrector[:,1]
            F_corrector[:,-1] = F_corrector[:,-2]
            S_corrector[:,0] = S_corrector[:,1]
            S_corrector[:,-1] = S_corrector[:,-2]

            Q = Q_corrector
            F = F_corrector
            S = S_corrector
            L2 = np.log(np.sqrt(np.average((rho_corrector - self.rho)**2)))
            L2s.append(L2)
            
            self.rho = rho_corrector
            self.u = u_corrector
            self.p = p_corrector
            self.e = e_corrector


            print("============================== Pressure ======================================")
            print (self.p)
            print("============================== Rho ===========================================")
            print (self.rho)
            print("============================== Speed =========================================")
            print (self.u)
            print("============================== Energy ========================================")
            print (self.e)
            # Mise à jour du plot en temps réel
            # Mise à jour des données pour la densité
            iteration += 1
            iterations.append(iteration)
            if iteration%10 == 0:
                line1.set_ydata(self.rho)
                ax1.relim()
                ax1.autoscale_view()

                # Mise à jour des données pour la pression
                line2.set_ydata(self.p)
                ax2.relim()
                ax2.autoscale_view()

                            
                # Mise à jour des données pour le nombre de mach
                line3.set_ydata(self.u/(self.gamma*self.p/self.rho)**0.5)
                ax3.relim()
                ax3.autoscale_view()

                line4.set_data(np.log(iterations), L2s)
                ax4.relim()
                ax4.autoscale_view()


                # Redessiner la figure

                plt.draw()
                plt.pause(0.01)

            live_time = live_time +  delta_t
            Delta_t = self.CFL/(abs(self.u) + (self.gamma*self.p/self.rho)**0.5)*self.mesh.delta_x
            delta_t = min(Delta_t)
            
            

        plt.ioff()
        plt.show()
    



    def Solve_Beam(self,theta = 0.5,ei = 0.3, ee = 0.1):
        def compute_flux(Q,xi):
            rhoA = Q[0]
            rhouA = Q[1]
            eA = Q[2]

            rho = rhoA/self.mesh.area[xi]
            u = rhouA / rhoA
            e = eA/self.mesh.area[xi]
            p = (self.gamma - 1) * (e - 0.5 * rho * u**2)
            
            H = (e + p) / rho  # Enthalpie totale

            # Flux vector F
            F = np.zeros_like(Q)
            F[0] = rhouA
            F[1] = (rho * u**2 + p)*self.mesh.area[xi]
            F[2] = u * (e + p) * self.mesh.area[xi]

            return F
        
        def compute_source(Q,xi):
            rhoA = Q[0]
            rhouA = Q[1]
            eA = Q[2]

            rho = rhoA/self.mesh.area[xi]
            u = rhouA / rhoA
            e = eA/self.mesh.area[xi]
            p = (self.gamma - 1) * (e - 0.5 * rho * u**2)
            
            H = (e + p) / rho  # Enthalpie totale

            # Flux vector F
            S = np.zeros_like(Q)
            S[0] = 0
            S[1] = p*self.mesh.dAdx[xi]
            S[2] = 0
            return S
        
        def compute_jacobian(Q,xi):
            rhoA = Q[0]
            rhouA = Q[1]
            eA = Q[2]

            rho = Q[0]/self.mesh.area[xi]
            u = rhouA / rhoA
            e = eA/self.mesh.area[xi]
            p = (self.gamma - 1) * (e - 0.5 * rho * u**2)
            H = (e + p) / rho  # Enthalpie totale

            # Calcul de la matrice Jacobienne A
            A = np.zeros((3, 3))
            A[0, 0] = 0
            A[0, 1] = 1
            A[0, 2] = 0

            A[1, 0] = (self.gamma - 3) / 2 * u**2
            A[1, 1] = (3 - self.gamma) * u
            A[1, 2] = self.gamma - 1

            A[2, 0] = u * ((self.gamma - 1) / 2 * u**2 - H)
            A[2, 1] = H - (self.gamma - 1) * u**2
            A[2, 2] = self.gamma * u

            return A
        def apply_boundary_conditions(Q):
            # For example, fixed boundary conditions tube !!!
            Q[:, 0] = Q[:, 2]      # Left boundary
            Q[:, 1] = Q[:, 2]
            Q[:, -2] = Q[:, -3]    # Right boundary
            Q[:, -1] = Q[:, -3]
            return Q
        def beam_warming_step(Q):
            dx = self.mesh.delta_x
            nx = self.mesh.n + 4  # Including ghost cells
            Q = apply_boundary_conditions(Q)
            Q_new = Q.copy()
            theta = 1  # Implicit weighting parameter

            # Dissipation coefficients
            epsilon_i = 0.25  # First-order dissipation coefficient
            epsilon_e = 0.01  # Fourth-order dissipation coefficient

            # Calculate primitive variables and maximum wave speed for CFL condition
            rhoA = Q[0, 2:-2]
            rhouA = Q[1, 2:-2]
            eA = Q[2, 2:-2]
            A = self.mesh.area
            rho = rhoA / A
            u = rhouA / rhoA
            e = eA / A
            p = (self.gamma - 1) * (e - 0.5 * rho * u**2)
            a = np.sqrt(self.gamma * p / rho)
            max_speed = np.max(np.abs(u) + a)
            dt = self.CFL * dx / max_speed

            for i in range(2, nx - 2):
                xi = i - 2  # Index for the area array without ghost cells

                # Compute fluxes at current and previous cells
                F_i = compute_flux(Q[:, i], xi)
                F_im1 = compute_flux(Q[:, i - 1], xi - 1)
                A_i = compute_jacobian(Q[:, i], xi)
                S_i = compute_source(Q[:, i], xi)
                # First-order dissipation term
                grad_Delta_x = (Q[:, i + 1] - Q[:, i - 1]) / (2 * dx)
                dissipation_term_1 = epsilon_i * grad_Delta_x

                # Fourth-order dissipation term
                grad2_Delta_x = (Q[:, i + 1] - 2 * Q[:, i] + Q[:, i - 1]) / (dx ** 2)
                dissipation_term_2 = epsilon_e * grad2_Delta_x

                # Right-Hand Side
                RHS = - (dt / dx) * (F_i - F_im1) - dt * S_i + dissipation_term_2 - dissipation_term_1

                # Left-Hand Side matrix including first-order dissipation
                LHS = np.eye(3) + (theta * dt / dx) * A_i 

                # Solve for dQ
                dQ = np.linalg.solve(LHS, RHS)

                # Update Q
                Q_new[:, i] += dQ

            return Q_new, dt
        


        plt.ion()
        fig, (ax1, ax2, ax3,ax4) = plt.subplots(4, 1, figsize=(8, 8))
        fig.suptitle("distribution of density, pressure ,mach number and L2, along x")
        # Initialiser les lignes pour la densité et la pression
        line1, = ax1.plot(self.mesh.x, self.rho, label='Densité ')
        ax1.set_xlabel('x')
        ax1.set_ylabel('Densité')
        ax1.legend()
        ax1.grid(True)

        line2, = ax2.plot(self.mesh.x, self.p, label='Pression ')
        ax2.set_xlabel('x')
        ax2.set_ylabel('Pression ')
        ax2.legend()
        ax2.grid(True)


        line3, = ax3.plot(self.mesh.x, self.u/(self.gamma*self.p/self.rho)**0.5, label='Mac Number ')
        ax3.set_xlabel('x')
        ax3.set_ylabel('Mach number')
        ax3.legend()
        ax3.grid(True)

        line4, = ax4.plot([], [], label='L2 Error')
        ax4.set_xlabel('log Iterations')
        ax4.set_ylabel('log L2 Error')
        ax4.legend()
        ax4.grid(True)


        live_time = 0
        iteration = 0 
         # Initial Q
        Q = np.zeros((3, self.mesh.n + 4))
        Q[0, 2:-2] = self.rho*self.mesh.area
        Q[1, 2:-2] = Q[0, 2:-2] * self.u
        Q[2, 2:-2] = self.e*self.mesh.area
        Q = apply_boundary_conditions(Q)

        live_time = 0
        while live_time < self.mesh.t:
            Q, dt = beam_warming_step(Q)
            # Apply boundary conditions
            Q = apply_boundary_conditions(Q)

            # Update live_time
            live_time += dt
            self.rho = Q[0,2:-2]/self.mesh.area
            self.u = Q[1,2:-2]/Q[0,2:-2]
            self.e = Q[2,2:-2]/self.mesh.area
            self.p = (self.gamma - 1) * (self.e - self.rho*self.u**2 * 0.5)
            print("============================== Pressure ======================================")
            print (self.p)
            print("============================== Rho ===========================================")
            print (self.rho)
            print("============================== Speed =========================================")
            print (self.u)
            print("============================== Energy ========================================")
            print (self.e)
            iterations = self.mesh.x
            if iteration%10 == 0:
                line1.set_ydata(self.rho)
                ax1.relim()
                ax1.autoscale_view()

                # Mise à jour des données pour la pression
                line2.set_ydata(self.p)
                ax2.relim()
                ax2.autoscale_view()

                            
                # Mise à jour des données pour le nombre de mach
                line3.set_ydata(self.u/(self.gamma*self.p/self.rho)**0.5)
                ax3.relim()
                ax3.autoscale_view()

                line4.set_data(np.log(iterations), (self.gamma*self.p/self.rho)**0.5)
                ax4.relim()
                ax4.autoscale_view()


                # Redessiner la figure

                plt.draw()
                plt.pause(0.01)
        
        plt.ioff()
        plt.show()


        

            


        





            
            


        


    



        
        

            
            
