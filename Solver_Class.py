import numpy as np
import Mesh_Class as Mesh
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

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
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 6))
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
            L2 = (np.average(rho_corrector)**2 - np.average(self.rho)**2)**0.5
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



                # Redessiner la figure

                plt.draw()
                plt.pause(0.01)

            live_time = live_time +  delta_t
            Delta_t = self.CFL/(abs(self.u) + (self.gamma*self.p/self.rho)**0.5)*self.mesh.delta_x
            delta_t = min(Delta_t)
            
            

        plt.ioff()
        plt.show()
    
    def Solve_Beam(self,theta = 0.5,ei = 0.3, ee = 0.1):
        #ititialisation
        I = np.identity(self.mesh.n + 2)
        Delta_t = self.CFL/(abs(self.u) + (self.gamma*self.p/self.rho)**0.5)*self.mesh.delta_x
        delta_t = min(Delta_t)
        A_matrix = []
        Q_matrix = []
        F_matrix = []
        S_matrix = []
        #initialisation de A avec ghost cells

        
        for i in range(0,self.mesh.n):
            H = (self.e[i] + self.p[i])/self.rho[i]
            A_matrix.append(np.array([[0,1,0],
                               [self.u[i]*((self.gamma-3)/2),(self.gamma-3)*self.u[i],self.gamma-1],
                               [self.u[i]*(((self.gamma-1)/2)*self.u[i]**2 - H),  H - (self.gamma-1)*self.u[i]**2, self.gamma*self.u[i]]]))
            
            Q_matrix.append(np.array([[self.rho[i]*self.mesh.area[i]
                             ,self.rho[i]*self.mesh.area[i]*self.u[i]
                             ,self.e[i]*self.mesh.area[i]]]))
            
            F_matrix.append(np.array([[self.rho[i]*self.u[i]**2*self.mesh.area[i],
                                       (self.rho[i]*self.u[i]**2 + self.p[i])*self.mesh.area[i],
                                       self.u[i]*self.mesh.area[i]*(self.e[i]*self.p[i])]]))
            
            S_matrix.append(np.array([[0,self.mesh.area[i]*self.mesh.dAdx,0]],dtype=object))
            
        A_left = A_matrix[0]
        A_right = A_matrix[-1]
        F_left = F_matrix[0]
        F_right = F_matrix[-1]
        S_left = S_matrix[0]
        S_right = S_matrix[-1]
        Q_left = Q_matrix[0]
        Q_right = Q_matrix[-1]




        A_matrix.append(A_right)
        A_matrix.insert(0,A_left)
        F_matrix.append(F_right)
        F_matrix.insert(0,F_left)
        S_matrix.append(S_right)
        S_matrix.insert(0,S_left)
        Q_matrix.append(Q_right)
        Q_matrix.insert(0,Q_left)

        





            
            


        


    



        
        

            
            
