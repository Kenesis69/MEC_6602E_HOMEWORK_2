
#include <iostream>
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


//RAW division OF 2 VECTORS --> VECTOR
std::vector<double> Vector_Raw_Divide(std::vector<double> A, std::vector<double> B){
    if (A.size() != B.size()) {cout << "Vector_Raw_Multiply dimension conflict" << endl; exit(1);}

    std::vector<double> C;

    for (int i = 0; i<A.size(); i++){
        C.push_back(A[i]/B[i]);
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

// substract vector A from B that fucking simple. is not commutatif
std::vector<double> Vector_Substract(std::vector<double> A,std::vector<double> B){
    if (A.size() != B.size()) {cout << "Vector_Raw_Multiply dimension conflict" << endl; exit(1);}

    std::vector<double> C;

    for (int i = 0; i< A.size();i++){
        C.push_back(A[i]-B[i]);
    }
    return C;
}

//addition vector A from B that fucking simple.
std::vector<double> Vector_Addition(std::vector<double> A,std::vector<double> B){
    if (A.size() != B.size()) {cout << "Vector_Raw_Multiply dimension conflict" << endl; exit(1);}

    std::vector<double> C;

    for (int i = 0; i< A.size();i++){
        C.push_back(A[i]+B[i]);
    }
    return C;
}



























// The class Mesh will have all the information regarding space and time parameters
class Mesh{
public:
    //variable
    double x1,x2,Delta_t,t,Delta_x;
    
    int n;
    std::vector<double> area; 

    //constructor
    Mesh(double a, double b, int c, int d, std::vector<double> e,double f ) : x1(a), x2(b), n(c), Delta_t(d), area(e), t(f) {
    Delta_x = (x2 - x1)/(n-1);
    cout <<"********************** Mesh *****************************" << std::endl;
    cout <<"beginning x = "<< Mesh::x1 << std::endl;
    cout <<"end x = "<< Mesh::x2 << std::endl;
    cout <<"Space step delta x = "<< Mesh::Delta_x << std::endl;
    cout <<"Time Step Delta_t = "<< Mesh::Delta_t << std::endl;
    cout <<"Time t = "<< Mesh::t << std::endl;
    cout <<"number of nodes n = "<< Mesh::n << std::endl;
    }; 


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
    double gamma = 1.40;
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
    };


    void Initialise_Euler_Element(){
        cout<< "INITIALATION MAIN VECTORS"<< endl;

        //------------------------INITIALISATION------------------------------//
        //Itermediate variable
        std::vector<double> rho_u = Vector_Raw_Multiply(init_rho,init_u); //simple ru
        std::vector<double> rho_u_squared = Vector_Raw_Multiply(rho_u,init_u); // simple ru squared
        

        // Initiasing dA over dx
        std::vector<double> dAdx;
        for (int i = 0; i< Mesh_to_use.n; i++){
            if (i == 0 || i == Mesh_to_use.n-1){
                dAdx.push_back(0);
            }
            else{
                dAdx.push_back((Mesh_to_use.area[i+1] - Mesh_to_use.area[i-1])/2);
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

        //initialisation, pas conformable utilise Q,F,S

        //Adding Ghost Cells (Fonctionne)
        cout <<"Cell Ghosting"<< endl;
        addColumns(Q,2,2);
        addColumns(F,2,2);
        addColumns(S,2,2);
        
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

        

        
        cout<< "Beginning of loop"<< endl;
        cout<< "Time Step "<< Mesh_to_use.Delta_t << endl;
        cout<< "Time "<< Mesh_to_use.t << endl;
        for (double k=0; k<Mesh_to_use.t; k= k + Mesh_to_use.Delta_t){
            cout<< "time: "<< k << endl;
            //predictor STEP
            for (int i =0; i <Q .size();i++){
                for (int j =2; j<Q[0].size()-2;j++){
                    //le gros du code est ici

                    //PREDICTOR STEP
                    Q_Predictor[i][j] = Q[i][j] - (Mesh_to_use.Delta_t / Mesh_to_use.Delta_x) * (F[i][j] - F[i][j - 1]) + Mesh_to_use.Delta_t * S[i][j];
                    
                    
                    
                    



                };
            };

            
            // Calcul 1 des propriétés
            for (int p = 0; p < Mesh_to_use.n; p++) {
                rho[p] = Q_Predictor[0][p + 2] / Mesh_to_use.area[p];
                u[p] = Q_Predictor[1][p + 2] / Q_Predictor[0][p + 2];
                e[p] = Q_Predictor[2][p + 2] / Mesh_to_use.area[p];
                P[p] = (gamma - 1) * (e[p] - 0.5 * rho[p] * u[p] * u[p]);
            }

 

            // calcul variation aire 

            std::vector<double> dAdx;
            for (int i = 0; i< Mesh_to_use.n; i++){
                if (i == 0 || i == Mesh_to_use.n-1){
                dAdx.push_back(0);
            }
                else{
                dAdx.push_back((Mesh_to_use.area[i-1] + Mesh_to_use.area[i+1])/2);
            }
        }


            // calcul de F et S
            F_Predictor[0] = Q_Predictor[1];
            F_Predictor[1] = Vector_Raw_Multiply(Vector_Addition(Vector_Raw_Multiply(rho,Vector_Raw_Square(u)),P),Mesh_to_use.area);
            F_Predictor[2] = Vector_Raw_Multiply(Vector_Raw_Multiply(Vector_Addition(e,P),u),Mesh_to_use.area);
            S_Predictor[0] = zero(Mesh_to_use.n);
            S_Predictor[1] = Vector_Raw_Multiply(dAdx,P);
            S_Predictor[2] = zero(Mesh_to_use.n);
            
        
            //PREDICTOR STEP
            for (int ii =0; ii <Q .size();ii++){
                for (int jj =2; jj<Q[0].size()-2;jj++){
                    

                    //PREDICTOR STEP
                    Q_Corrected[ii][jj] = Q[ii][jj] - Mesh_to_use.Delta_t*(F_Predictor[ii][jj+1] - F_Predictor[ii][jj])/Mesh_to_use.Delta_x + Mesh_to_use.Delta_t*S_Predictor[ii][jj];
                    

                    
                    




                };
            };


            // Calcul 2 des propriétés
            for (int p = 0; p < Mesh_to_use.n; p++) {
                rho[p] = Q_Corrected[0][p + 2] / Mesh_to_use.area[p];
                u[p] = Q_Corrected[1][p + 2] / Q_Corrected[0][p + 2];
                e[p] = Q_Corrected[2][p + 2] / Mesh_to_use.area[p];
                P[p] = (gamma - 1) * (e[p] - 0.5 * rho[p] * u[p] * u[p]);
            }

            F_Corrected[0] = Q_Predictor[1];
            F_Corrected[1] = Vector_Raw_Multiply(Vector_Addition(Vector_Raw_Multiply(rho,Vector_Raw_Square(u)),P),Mesh_to_use.area);
            F_Corrected[2] = Vector_Raw_Multiply(Vector_Raw_Multiply(Vector_Addition(e,P),u),Mesh_to_use.area);
            S_Corrected[0] = zero(Mesh_to_use.n);
            S_Corrected[1] = Vector_Raw_Multiply(dAdx,P);
            S_Corrected[2] = zero(Mesh_to_use.n);


            F = F_Corrected;
            S = S_Corrected;
            Q = Q_Corrected;
            printMatrix(Q_Corrected);

            





           

            
        };

        
    };



    
};














//Main
int main(){
    //definition parameters
    double c = 1; //constant c
    double x1 = 0.0; // beginning
    double x2  = 1000 ; // end
    int n = 4; // number of point
    std::vector<double> area; // definition de l'aire
    string method = "" ;

    double delta_x = (x2 - x1)/(n-1); //delta_x calculated
    double delta_t = 1; //delta_t 
    double t = 250;
    double CFL = (c*delta_t)/delta_x; // CFL calculated
    cout <<"********************** PARAMETERS ***************************" << std::endl;
    cout<< "The CFL is:  "<<CFL << std::endl; 

    //defining area for now is case 1
    for (int i = 0; i<n ; i++){
        area.push_back(1.0);
    }    



    // Beginning the Solve Process
    Mesh Mesh_devoir(x1,x2,n,delta_t,area,t); // declaration object Mesh
    
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
            Init_e.push_back(0);
        
        } else {
            Init_p.push_back(4);
            Init_rho.push_back(4);
            Init_u.push_back(0);
            Init_e.push_back(0);
            
        }
    }
    Solveur_Devoir.Condition_initiales(Init_p,Init_rho,Init_u,Init_e);
    Solveur_Devoir.Initialise_Euler_Element();
    Solveur_Devoir.MacCormack();
    

    return 0;
    
};