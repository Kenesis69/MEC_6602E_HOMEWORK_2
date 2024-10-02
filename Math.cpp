#include <iostream>
#include <vector>
using namespace std; 
#include "Math.h"

class math
{
private:
    /* data */
public:
   
  
    std::vector<double> A;
    std::vector<double> B;

    std::vector<double> Simple_Multiply(){
        if (A.size() != B.size()){
            cout<< "Vecteur pas de la bonne taille"<< endl;
        }
        std::vector<double> C;
        for (int i = 0; i<A.size(); i++){
            C.push_back(A[i]*B[i]);
        }
        return C;
    }
    
};



