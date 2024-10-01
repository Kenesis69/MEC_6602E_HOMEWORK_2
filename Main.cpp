
#include <iostream>
#include <vector>
using namespace std; 
constexpr double pi = 3.14159265358979323846;

// The class Mesh will have all the information regarding space and time parameters
class Mesh{
public:
    //variable
    double x1,x2,Delta_t;
    int n;
    

    //constructor
    Mesh(double a, double b, int c, int d) : x1(a), x2(b), n(c), Delta_t(d) {
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
        cout << "Vector is: " << std::endl;
        for (int i = 0; i < result.size(); i++) {
            cout << result.at(i) << " ";  
        }
        cout << endl;
        return result;
    };
};

class Solveur{
public:
    Mesh Mesh_to_use;
    string method;
    std::vector<double> init_p, init_rho, init_u;

    // Constructor
    Solveur (Mesh a, string b) : Mesh_to_use(a), method(b) {
    cout <<"********************** SOLVEUR *****************************" << std::endl;
    cout << "Mesh used :" << endl;
    };

    // Creation matrice 1D contenant les conditions initiales 
    void Condition_initiales(std::vector<double> p, std::vector<double>rho,std::vector<double>u){
                if (p.size() != Mesh_to_use.n || rho.size() != Mesh_to_use.n || u.size() != Mesh_to_use.n){
                    cout<< "ERROR: VECTOR DOES NOT MATCH NUMBHER OF ELEMENT"<< endl;
                    exit(1);
                }
                else {cout<< "Condition has been succesfully entered"<< endl;}
                init_p = p;
                init_rho = rho;
                init_u = u;
    };

    void Solve_Mac_Cormack_Method(){
        cout<< "Beginning Mac Cormack Method"<< endl;
        // initialisation des matrice Q,E,S selon les conditions initiales
        
    };

    
};

//Main
int main(){
    //definition parameters
    double c = 0.5; //constant c
    double x1 = 0.0; // beginning
    double x2  = 1000 ; // end
    int n = 5; // number of point
    double A = 1.0; // definition de l'aire
    string method = "" ;

    double delta_x = (x2 - x1)/(n-1); //delta_x calculated
    double delta_t = 2; //delta_t 
    double CFL = (c*delta_t)/delta_x; // CFL calculated
    cout <<"********************** PARAMETERS ***************************" << std::endl;
    cout<< "The CFL is:  "<<CFL << std::endl; 


    // Beginning the Solve Process
    Mesh Mesh_devoir(x1,x2,n,delta_t); // declaration object Mesh
    
    Solveur Solveur_Devoir(Mesh_devoir,method); // Declaration solveur
    
    // Initial conditions
    std::vector<double> Init_p ;
    std::vector<double> Init_rho ;
    std::vector<double> Init_u ;
    std::vector<double> x_vector = Mesh_devoir.get_vector();
    for (int i = 0; i < 5;i++) {
        if (x_vector[i] < 500 ){
            Init_p.push_back(1);
            Init_rho.push_back(1);
            Init_u.push_back(0);
        
        } else {
            Init_p.push_back(4);
            Init_rho.push_back(4);
            Init_u.push_back(0);
        }
    }

    Solveur_Devoir.Condition_initiales(Init_p,Init_rho,Init_u);



    return 0;
    
};