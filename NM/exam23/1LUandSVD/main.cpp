#include <iostream>
#include <fstream>
#include "nr3.h"
#include "utils.h"
#include "cholesky.h"
#include "ludcmp.h"
#include <cmath>
#include "svd.h"

//your extra include headers

using namespace std;

void print_error(SVD SVD_A){
    MatDoub V = SVD_A.v;
    VecDoub W = SVD_A.w;

    VecDoub dx(W.size());
    double sum;
    for (int i = 0; i < W.size(); i++)
    {
        sum = 0;
        for (int j = 0; j < W.size(); j++)
        {
            if (W[i] > pow(10, -15)){
                sum += pow(V[j][i]/W[i], 2);
            }
        }
        dx[i] = sqrt(sum);
    }
    utils::print(dx, "delta x");
}


int main() {

    ifstream input_A("Ex1A.dat");
    int Ma, Na;
    input_A >> Ma;
    input_A >> Na;
    cout << Ma << " " << Na << endl;
    MatDoub A(Ma,Na);
    for (int i = 0; i < Ma; i++)
    {
        for (int j = 0; j < Na; j++)
        {
            input_A >> A[i][j];
        }
    }
    

    ifstream input_b("Ex1b.dat");
    int Mb, Nb;
    input_b >> Mb;
    input_b >> Nb;
    // cout << Mb << " " << Nb << endl;
    VecDoub b(Mb);
    for (int i = 0; i < Mb; i++)
    {
        input_b >> b[i];
    }    

    // 1.

    MatDoub C(Ma,Na);
    C = utils::Transpose(A)*A;
    VecDoub c(Na);
    c = utils::Transpose(A)*b;

    utils::printDiag(C, "C");
    utils::print(c, "c");

    // 2

    LUdcmp lu(C);

    // utils::print(lu.lu, "LU"); // diagonal og over = U, alt under diagonal = L
    utils::printDiag(lu.lu, "U diagonal elements");

    cout << "pivoting bookkeeping. Ikke helt sikker pÅ hvordan det skal tolkes" << endl;
    for (int i = 0; i < Na; i++)
    {
        cout << lu.indx[i] << " ";
    }
    cout << endl;

    // 3
    VecDoub x(Na);
    lu.solve(c, x);
    utils::print(x, "x");


    // 4
    /*
the componentwise error in the solution xx can be estimated using Eq. 15.4.19 from the SVD formulation:
    Eq. 15.4.19 
    Fra lek 4 slide 18
This expression provides insight into how the condition of each singular direction affects the accuracy of each component of xx, independently of the right-hand side bb.
    
    */

    return 0;

}
