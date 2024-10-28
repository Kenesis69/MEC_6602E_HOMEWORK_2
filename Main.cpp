#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
using namespace std; 
constexpr double pi = 3.14159265358979323846;


//print fonction
template <typename T>
void printMatrix(const vector<vector<T>>& matrix) {
    for (const auto& row : matrix) {
        for (const T& element : row) {
            cout << element << " ";
        }
        cout << endl; // Move to the next line after printing each row
    }
    cout<< "----------------------------------"<< endl;
}

// Function to print a vector
template <typename T>
void printVector(const std::vector<T>& vec) {
    for (const T& element : vec) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
}



void removeColumns(std::vector<std::vector<double>>& matrix, int numColumnsBegin, int numColumnsEnd) {
    // Supprimer des colonnes au début et à la fin pour chaque ligne
    for (auto& row : matrix) {
        // Vérifier que la taille de la ligne est suffisamment grande
        if (row.size() > numColumnsBegin + numColumnsEnd) {
            // Supprimer les colonnes du début
            row.erase(row.begin(), row.begin() + numColumnsBegin);

            // Supprimer les colonnes de la fin
            row.erase(row.end() - numColumnsEnd, row.end());
        } else {
            std::cerr << "Error: Not enough columns to remove" << std::endl;
        }
    }
}

void addColumns(std::vector<std::vector<double>>& matrix, int numColumnsBegin, int numColumnsEnd) {
    // Add columns at the beginning and the end for each row
    for (auto& row : matrix) {
        // Add columns at the beginning (shift elements to the right)
        row.insert(row.begin(), numColumnsBegin, row.front());

        // Add columns at the end
        row.insert(row.end(), numColumnsEnd, row.back());
    }
}












//MATH FONCTIONS
//zeros fonction np.zero
std::vector<double> zero(int n){
    std::vector<double> C;
    for (int i = 0; i < n; i++){
        C.push_back(0);
    }
    return C;
}

//RAW MULTIPLICATION OF 2 VECTORS --> VECTOR
std::vector<double> Vector_Raw_Multiply(std::vector<double> A, std::vector<double> B){
    if (A.size() != B.size()) {cout << "Vector_Raw_Multiply dimension conflict" << endl; exit(1);}

    std::vector<double> C;

    for (int i = 0; i<A.size(); i++){
        C.push_back(A[i]*B[i]);
    };
    return C;
};


//RAW DIVISION OF 2 VECTORS --> VECTOR STUPID WILL INVERSE THE ORDER
std::vector<double> Vector_Raw_Divide(std::vector<double> A, std::vector<double> B){
    if (A.size() != B.size()) {cout << "Vector_Raw_Multiply dimension conflict" << endl; exit(1);}

    std::vector<double> C;

    for (int i = 0; i<A.size(); i++){
        C.push_back(A[i]/B[i]);
    };
    return C;
};

// Multiply Scalar of a vector
std::vector<double> Vector_Scalar_Multiply(double A, std::vector<double> B){
    
    std::vector<double> C;

    for (int i = 0; i<B.size(); i++){
        C.push_back(A*B[i]);
    };
    return C;
};

//vector of scalar Divided by vector

std::vector<double> Vector_Scalar_Divided_by_vector(double A, std::vector<double> B){

    std::vector<double> C;

    for (int i = 0; i<B.size(); i++){
        C.push_back(A/B[i]);
    };
    return C;
};

//RAW SQUARIFICATION OF 1 VECTORS --> VECTOR
std::vector<double> Vector_Raw_Square(std::vector<double> A){
    

    std::vector<double> C;

    for (int i = 0; i<A.size(); i++){
        C.push_back(A[i]*A[i]);
    };
    return C; 
};

std::vector<double> Vector_Raw_racine(std::vector<double> A){
    

    std::vector<double> C;

    for (int i = 0; i<A.size(); i++){
        C.push_back(std::sqrt(A[i]));
    };
    return C; 
};

std::vector<double> Vector_Raw_Abs(std::vector<double> A){
    

    std::vector<double> C;

    for (int i = 0; i<A.size(); i++){
        C.push_back(std::abs(A[i]));
    };
    return C; 
};


//addition vector A from B that fucking simple.
std::vector<double> Vector_Addition(std::vector<double> A,std::vector<double> B){
    if (A.size() != B.size()) {cout << "Vector_Raw_Multiply dimension conflict" << endl; exit(1);}

    std::vector<double> C;

    for (int i = 0; i< A.size();i++){
        C.push_back(A[i]+B[i]);
    }
    return C;
}



double Vector_Raw_Min(std::vector<double> A){
    

    double C = A[0];

    for (int i = 1; i<A.size(); i++){
        if (A[i] < C){
            C = A[i];
        };
    };
    return C; 
};

























// The class Mesh will have all the information regarding space and time parameters
class Mesh{
public:
    //variable
    double x1,x2,Delta_t,t,Delta_x;
    
    int n;
    std::vector<double> area; 
    std::vector<double> dAdx;
    //constructor
    Mesh(double a, double b, int c, std::vector<double> e,double f ) : x1(a), x2(b), n(c), area(e), t(f) {
    Delta_x = (x2 - x1)/(n-1);
    cout <<"********************** Mesh *****************************" << std::endl;
    cout <<"beginning x = "<< Mesh::x1 << std::endl;
    cout <<"end x = "<< Mesh::x2 << std::endl;
    cout <<"Space step delta x = "<< Mesh::Delta_x << std::endl;
    cout <<"Time Step Delta_t = "<< Mesh::Delta_t << std::endl;
    cout <<"Time t = "<< Mesh::t << std::endl;
    cout <<"number of nodes n = "<< Mesh::n << std::endl;
    cout <<"AREA IS "<< std::endl;
    printVector(area);
    
    for (int i = 0; i < n; i++){
        if (i == 0){
            // Forward difference at the left boundary
            dAdx.push_back((area[1] - area[0]) / Delta_x);
        }
        else if (i == n - 1){
            // Backward difference at the right boundary
            dAdx.push_back((area[n - 1] - area[n - 2]) / Delta_x);
        }
        else{
            // Central difference
            dAdx.push_back((area[i + 1] - area[i - 1]) / (2 * Delta_x));
        }
    }
    }; 

                // calcul variation aire 

            



    // vector getter and printer
    std::vector<double> get_vector(){
    vector<double> result;
    
        for(int i = 0; i < n ; i++){
            result.push_back(x1 + (i)*Delta_x);
        };
        //cout << "Vector is: " << std::endl;
        //for (int i = 0; i < result.size(); i++) {
            //cout << result.at(i) << " ";  
        //}
        cout << endl;
        return result;
    };
};































class Solveur{
public:

    
    Mesh Mesh_to_use;
    string method;
    std::vector<double> init_p, init_rho, init_u, init_e;
    double gamma = 1.4;
    std::vector<std::vector<double>> Q, F, S;
     

    // Constructor
    Solveur (Mesh a, string b) : Mesh_to_use(a), method(b) {
    cout <<"********************** SOLVEUR *****************************" << std::endl;
    };

    // Creation matrice 1D contenant les conditions initiales 
    void Condition_initiales(std::vector<double> p, std::vector<double>rho,std::vector<double>u, std::vector<double> e){
                if (p.size() != Mesh_to_use.n || rho.size() != Mesh_to_use.n || u.size() != Mesh_to_use.n){
                    cout<< "ERROR: VECTOR DOES NOT MATCH NUMBHER OF ELEMENT"<< endl;
                    exit(1);
                }
                else {cout<< "Condition has been succesfully entered"<< endl;}
                init_p = p;
                init_rho = rho;
                init_u = u;
                init_e = e;
                cout<< "Here are the initials conditions"<< endl;
                std::vector<std::vector<double>> init_cond = {init_rho,init_u,init_p};
                printMatrix(init_cond);
    };


    void Initialise_Euler_Element(){
        cout<< "INITIALATION MAIN VECTORS"<< endl;

        //------------------------INITIALISATION------------------------------//
        //Itermediate variable
        std::vector<double> rho_u = Vector_Raw_Multiply(init_rho,init_u); //simple ru
        std::vector<double> rho_u_squared =Vector_Raw_Multiply(init_rho,Vector_Raw_Square(init_u)) ; // simple ru squared
        

        // Initiasing dA over dx
        
        std::vector<double> dAdx;
        for (int i = 0; i < Mesh_to_use.n; i++){
            if (i == 0){
                // Forward difference at the left boundary
                dAdx.push_back((Mesh_to_use.area[1] - Mesh_to_use.area[0]) / Mesh_to_use.Delta_x);
            }
            else if (i == Mesh_to_use.n - 1){
                // Backward difference at the right boundary
                dAdx.push_back((Mesh_to_use.area[Mesh_to_use.n - 1] - Mesh_to_use.area[Mesh_to_use.n - 2]) / Mesh_to_use.Delta_x);
            }
            else{
                // Central difference
                dAdx.push_back((Mesh_to_use.area[i + 1] - Mesh_to_use.area[i - 1]) / (2 * Mesh_to_use.Delta_x));
            }
        }


        // initialisation des matrice Q,E,S selon les conditions initiales
        Q = {Vector_Raw_Multiply(init_rho,Mesh_to_use.area),Vector_Raw_Multiply(rho_u,Mesh_to_use.area),Vector_Raw_Multiply(init_e,Mesh_to_use.area)};
        F = {Vector_Raw_Multiply(rho_u,Mesh_to_use.area),Vector_Raw_Multiply(Mesh_to_use.area,Vector_Addition(rho_u_squared,init_p)),Vector_Raw_Multiply(Vector_Raw_Multiply(Mesh_to_use.area,Vector_Addition(init_e,init_p)),init_u)};
        S = {zero(Mesh_to_use.n),Vector_Raw_Multiply(init_p,dAdx),zero(Mesh_to_use.n)};
        cout<<"Q initial"<< endl;
        printMatrix(Q);
        cout<<"F initial"<< endl;
        printMatrix(F);
        cout<<"S initial"<< endl;
        printMatrix(S);
        
    };

    void MacCormack(){

        

        //Adding Ghost Cells (Fonctionne)
        cout <<"Cell Ghosting"<< endl;
        addColumns(Q,1,1);
        addColumns(F,1,1);
        addColumns(S,1,1);
        
        printMatrix(Q);
        printMatrix(F);
        printMatrix(S);

        //Début de la méthode MacCormack
        double live_time = 0;
        std::vector<std::vector<double>> Q_Corrected = Q;
        std::vector<std::vector<double>> Q_Predictor = Q;
        std::vector<std::vector<double>> F_Corrected = F;
        std::vector<std::vector<double>> S_Corrected = S;
        std::vector<std::vector<double>> F_Predictor = F;
        std::vector<std::vector<double>> S_Predictor = S;

        std::vector<double> rho(Mesh_to_use.n);
        std::vector<double> P(Mesh_to_use.n);
        std::vector<double> e(Mesh_to_use.n);
        std::vector<double> u(Mesh_to_use.n);

        std::vector<double> rho_Predictor(Mesh_to_use.n);
        std::vector<double> P_Predictor(Mesh_to_use.n);
        std::vector<double> e_Predictor(Mesh_to_use.n);
        std::vector<double> u_Predictor(Mesh_to_use.n);

        std::vector<double> rho_Corrected(Mesh_to_use.n);
        std::vector<double> P_Corrected(Mesh_to_use.n);
        std::vector<double> e_Corrected(Mesh_to_use.n);
        std::vector<double> u_Corrected(Mesh_to_use.n);
 


        std::vector<std::vector<double>> result;
        std::vector<double> Delta_t;
        Delta_t = Vector_Scalar_Divided_by_vector(Mesh_to_use.Delta_x,Vector_Addition(Vector_Raw_Abs(init_u) ,Vector_Raw_racine(Vector_Scalar_Multiply(gamma, Vector_Raw_Divide(init_p,init_rho)))));;
        double delta_t = 0.5*Vector_Raw_Min(Delta_t);
        cout<< "Beginning of loop"<< endl;
        
        cout<< "Time "<< Mesh_to_use.t << endl;
        for (double k=0; k<Mesh_to_use.t; k= k + delta_t){
            cout<< "time: "<< k << endl;
            //predictor STEP
            for (int i =0; i <Q .size();i++){
                for (int j =1; j<Mesh_to_use.n +1;j++){
                   

                    
                    // Predictor Step (Corrected forward difference)
                    Q_Predictor[i][j] = Q[i][j] - ((delta_t) * (F[i][j+1] - F[i][j])/Mesh_to_use.Delta_x) + delta_t * S[i][j];


                    
                    
                    
                    



                };
            };



            // Calcul predictif des propriétés
            for (int p = 0; p < Mesh_to_use.n; p++) {
                rho_Predictor[p] = Q_Predictor[0][p + 1] / Mesh_to_use.area[p];
                u_Predictor[p] = Q_Predictor[1][p + 1] / Q_Predictor[0][p + 1];
                e_Predictor[p] = Q_Predictor[2][p + 1] / Mesh_to_use.area[p];
                P_Predictor[p] = (gamma - 1) * (e_Predictor[p] - 0.5 * rho_Predictor[p] * u_Predictor[p] * u_Predictor[p]);
                if (P_Predictor[p] < 0){
                    P_Predictor[p] = 0;
                }
            }

 




            // calcul de F et S
            F_Predictor[0] = Q_Predictor[1];
            F_Predictor[1] = Vector_Raw_Multiply(Vector_Addition(Vector_Raw_Multiply(rho_Predictor,Vector_Raw_Square(u_Predictor)),P_Predictor),Mesh_to_use.area);
            F_Predictor[2] = Vector_Raw_Multiply(Vector_Raw_Multiply(Vector_Addition(e_Predictor,P_Predictor),u_Predictor),Mesh_to_use.area);
            S_Predictor[0] = zero(Mesh_to_use.n);
            S_Predictor[1] = Vector_Raw_Multiply(Mesh_to_use.dAdx,P_Predictor);
            S_Predictor[2] = zero(Mesh_to_use.n);
            
                        int total_cells = Mesh_to_use.n + 2; // Total cells including ghost cells
            for (int var = 0; var < Q.size(); var++) {
                // Left boundary (indices 0 and 1)
                Q_Predictor[var][0] = Q_Predictor[var][1];
                F_Predictor[var][0] = F_Predictor[var][1];
                S_Predictor[var][0] = S_Predictor[var][1];
                // Right boundary (indices total_cells - 2 and total_cells - 1)
                Q_Predictor[var][total_cells - 1] = Q_Predictor[var][total_cells - 2];
                F_Predictor[var][total_cells - 1] = F_Predictor[var][total_cells - 2];
                S_Predictor[var][total_cells - 1] = S_Predictor[var][total_cells - 2];
            }

            //Correcteur STEP
            for (int k =0; k <Q .size();k++){
                for (int m =1; m < Mesh_to_use.n +1;m++){
                    

                  
                    // Corrector Step 
                    Q_Corrected[k][m] = Q[k][m]  - ((delta_t) * (F_Predictor[k][m] - F_Predictor[k][m-1]) / Mesh_to_use.Delta_x) + delta_t * S_Predictor[k][m];
                    
                };
            };

            for (int n =0; n <Q .size();n++){
                for (int o =1; o<Mesh_to_use.n +1;o++){
                    // updading
                    Q_Corrected[n][o] = 0.5*(Q_Predictor[n][o] + Q_Corrected[n][o]);
                };
            };





            // Apply boundary conditions to Q_Corrected
        
             // Total cells including ghost cells
            for (int var = 0; var < Q.size(); var++) {
                // Left boundary (indices 0 and 1)
                
                Q_Corrected[var][0] = Q_Corrected[var][1];
                F_Corrected[var][0] = F_Corrected[var][1];
                S_Corrected[var][0] = S_Corrected[var][1];
                // Right boundary (indices total_cells - 2 and total_cells - 1)
                Q_Corrected[var][total_cells - 1] = Q_Corrected[var][total_cells - 2];
                F_Corrected[var][total_cells - 1] = F_Corrected[var][total_cells - 2];
                S_Corrected[var][total_cells - 1] = S_Corrected[var][total_cells - 2];
            }
            
                        // Calcul correctif des proprietés
            for (int p = 0; p < Mesh_to_use.n; p++) {
                rho_Corrected[p] = Q_Corrected[0][p + 1] / Mesh_to_use.area[p];
                u_Corrected[p] = Q_Corrected[1][p + 1] / Q_Corrected[0][p + 1];
                e_Corrected[p] = Q_Corrected[2][p + 1] / Mesh_to_use.area[p];
                P_Corrected[p] = (gamma - 1) * (e_Corrected[p] - 0.5 * rho_Corrected[p] * u_Corrected[p] * u_Corrected[p]);
                if (P_Corrected[p] < 0){
                    P_Corrected[p] = 0;
                }
            }

            F_Corrected[0] = Q_Corrected[1];
            F_Corrected[1] = Vector_Raw_Multiply(Vector_Addition(Vector_Raw_Multiply(rho_Corrected,Vector_Raw_Square(u_Corrected)),P_Corrected),Mesh_to_use.area);
            F_Corrected[2] = Vector_Raw_Multiply(Vector_Raw_Multiply(Vector_Addition(e_Corrected,P_Corrected),u_Corrected),Mesh_to_use.area);
            S_Corrected[0] = zero(Mesh_to_use.n);
            S_Corrected[1] = Vector_Raw_Multiply(Mesh_to_use.dAdx,P_Corrected);
            S_Corrected[2] = zero(Mesh_to_use.n);
            
            Q = Q_Corrected;
            F = F_Corrected;
            S = S_Corrected;

            rho = rho_Corrected;
            u = u_Corrected;
            P = P_Corrected;
            e = e_Corrected;

            //delta_t/delta_x)*(abs(Init_u[i]) + std::sqrt(gamma*Init_p[i]/Init_rho[i] < 1
            Delta_t = Vector_Scalar_Divided_by_vector(Mesh_to_use.Delta_x,Vector_Addition(Vector_Raw_Abs(u) ,Vector_Raw_racine(Vector_Scalar_Multiply(gamma, Vector_Raw_Divide(P,rho)))));
            delta_t = Vector_Raw_Min(Delta_t)/2;
            cout << "delta_t min = "<< delta_t<< endl;
            result = {rho_Corrected, u_Corrected,e_Corrected, P_Corrected};
            printMatrix(result);
            
            
            
            //cout << "---------------------Q IS---------------------------------------------"<<endl;
            //printMatrix(Q);
            //cout << "---------------------F IS---------------------------------------------"<<endl;
            //printMatrix(F);
            //cout << "---------------------S IS---------------------------------------------"<<endl;
            //printMatrix(S);
        
           


            





           

            
        };

        
    };

    void Beam_Warming(){

    }



    
};














//Main
int main(){
    //definition parameters
    double gamma = 1.4;
    double x1 = 0.0; // beginning
    double x2  = 1000.0 ; // end
    int n = 21; // number of point
    std::vector<double> area; // definition de l'aire
    string method = "" ;

    double t = 250.0;



    //defining area for now in case 1
    for (int i = 0; i<n ; i++){
        area.push_back(1.0);
    }    



    // Beginning the Solve Process
    Mesh Mesh_devoir(x1,x2,n,area,t); // declaration object Mesh
    
    Solveur Solveur_Devoir(Mesh_devoir,method); // Declaration solveur
    
    // Initial conditions
    std::vector<double> Init_p;
    std::vector<double> Init_rho;
    std::vector<double> Init_u;
    std::vector<double> Init_e;
    std::vector<double> x_vector = Mesh_devoir.get_vector();

    
    for (int i = 0; i < n;i++) {
        if (x_vector[i] < 500 ){
            Init_p.push_back(1);
            Init_rho.push_back(1);
            Init_u.push_back(0);
            Init_e.push_back(Init_p[i]/((gamma - 1) ));
        
        } else {
            Init_p.push_back(0.125);
            Init_rho.push_back(0.1);
            Init_u.push_back(0);
            Init_e.push_back(Init_p[i] / ((gamma - 1)));
        }
    }
    

    Solveur_Devoir.Condition_initiales(Init_p,Init_rho,Init_u,Init_e);
    Solveur_Devoir.Initialise_Euler_Element();
    Solveur_Devoir.MacCormack();
    

    return 0;
    
};

