import numpy as np


class Mesh():
    def __init__(self, x1, x2, n, area, t):
        self.x1 = x1
        self.x2 = x2
        self.n = n
        self.area = area
        self.t = t
        #ajout delta_x
        self.x = np.linspace(x1,x2,n)
        delta_x = (x2 - x1)/(n-1)
        
        #Ajout variation d'aire
        dAdx = np.zeros(n)
        for i in range(0,self.n):
            dAdx[i] = 0.2776*np.cosh(4.0-0.8*self.x[i])**-2 #tuyere
            #dAdx[i] = 0 #tube
            #if i == 0:
                #dAdx[i] = (area[i+1]-area[i])/delta_x
            #elif i == n-1:
                #dAdx[i] = (area[i]-area[i-1])/delta_x
            #else : dAdx[i] = (area[i+1]-area[i-1])/(2*delta_x)

        
        self.dAdx = dAdx
        self.delta_x = delta_x
        print("********************** Maillage *****************************")
        print(f"beginning x = {self.x1}")
        print(f"end x = {self.x2}")
        print(f"Space step delta x = {self.delta_x}")
        print(f"Time t = {self.t}")
        print(f"number of nodes n = {self.n}")
        print("AREA IS")
        print(area)
        print("AREA VARIATION IS")
        print(dAdx)

        

        
        

            
                

