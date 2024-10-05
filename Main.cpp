
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
    double x1,x2,Delta_t,t;
    int n;
    std::vector<double> area; 

    //constructor
    Mesh(double a, double b, int c, int d, std::vector<double> e,double f ) : x1(a), x2(b), n(c), Delta_t(d), area(e), t(f) {
    cout <<"********************** Mesh *****************************" << std::endl;
    cout <<"beginning x = "<< Mesh::x1 << std::endl;
    cout <<"end x = "<< Mesh::x2 << std::endl;
    cout <<"Time Step Delta_t = "<< Mesh::Delta_t << std::endl;
    cout <<"number of nodes n = "<< Mesh::n << std::endl;
    }; 


    // vector getter and printer
    std::vector<double> get_vector(){
        vector<double> result;
        double delta_x = (x2 - x1)/(n-1);
        for(int i = 0; i < n ; i++){
            result.push_back(x1 + (i)*delta_x);
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
        

        // now we initialise P = (gamma - 1)(e - pu^2/2)
        std::vector<double> P;
        for (int i = 0; i< Mesh_to_use.n; i++){
            P.push_back((gamma-1)*(init_e[i]- rho_u_squared[i]*0.5));
        }
        

        // Initiasing dA over dx
        std::vector<double> dAdx;
        for (int i = 0; i< Mesh_to_use.n; i++){
            if (i == 0 || i == Mesh_to_use.n-1){
                dAdx.push_back(0);
            }
            else{
                dAdx.push_back((Mesh_to_use.area[i-1] + Mesh_to_use.area[i+1])/2);
            }
        }


        // initialisation des matrice Q,E,S selon les conditions initiales
        Q = {Vector_Raw_Multiply(init_rho,Mesh_to_use.area),Vector_Raw_Multiply(rho_u,Mesh_to_use.area),Vector_Raw_Multiply(init_e,Mesh_to_use.area)};
        F = {Vector_Raw_Multiply(rho_u,Mesh_to_use.area),Vector_Raw_Multiply(Mesh_to_use.area,Vector_Addition(rho_u_squared,P)),Vector_Raw_Multiply(Vector_Raw_Multiply(Mesh_to_use.area,Vector_Addition(init_e,P)),init_u)};
        S = {zero(Mesh_to_use.n),Vector_Raw_Multiply(P,dAdx),zero(Mesh_to_use.n)};
        cout<<"Q initial"<< endl;
        printMatrix(Q);
        
    };

    void MacCormack(){

        //initialisation, pas conformable utilise Q,F,S

        //Adding Ghost Cells (Fonctionne)
        cout <<"Cell Ghosting"<< endl;
        addColumns(Q,2,2);
        addColumns(F,2,2);
        addColumns(S,2,2);
        printMatrix(Q);

        //Début de la méthode MacCormack
        double live_time = 0;
        std::vector<std::vector<double>> Q_future = Q;
        std::vector<std::vector<double>> F_future = F;
        std::vector<std::vector<double>> S_future = S;

        

        while (live_time < Mesh_to_use.t){

            //Double Boucle bébé
            for (int i =0; i<Q.size();i++){
                for (int j =2; i<Q[0].size()-2;i++){
                    //le gros du code est ici

                    //PREDICTOR STEP
                    Q_future[i][j] = Q[i][j] - Mesh_to_use.Delta_t*(F[i][j]-F[i][j-1]) + Mesh_to_use.Delta_t*S[i][j];
                    printMatrix(Q_future);

                    //CORRECTORSTEP



















                }
            }
            live_time +=  Mesh_to_use.Delta_t;
        };

        
    }



    
};














//Main
int main(){
    //definition parameters
    double c = 0.5; //constant c
    double x1 = 0.0; // beginning
    double x2  = 1000 ; // end
    int n = 5; // number of point
    std::vector<double> area; // definition de l'aire
    string method = "" ;

    double delta_x = (x2 - x1)/(n-1); //delta_x calculated
    double delta_t = 1; //delta_t 
    double t = 1000;
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