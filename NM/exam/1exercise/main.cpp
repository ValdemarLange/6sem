#include <iostream>
#include <fstream>
#include "nr3.h"
#include "utils.h"
#include <cmath>
#include "svd.h"

using namespace std;

int main() {

    ifstream input_A("NUM S25_Ex1A.dat");
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
    

    ifstream input_b("NUM S25_Ex1b.dat");
    int Mb, Nb;
    input_b >> Mb;
    input_b >> Nb;
    cout << Mb << " " << Nb << endl;
    VecDoub b(Mb);
    for (int i = 0; i < Mb; i++)
    {
        input_b >> b[i];
    }    

    // 1. Find the Singular Value Decomposition A = UWV^T . State the diagonal elements in W
    SVD svd(A);
    utils::print(svd.w, "w diagonal elements");

    // 2. here is a single element in W that is basically zero. Use the information
    // from the SVD matrices to state a unit vector in the null space of A.
    VecDoub null_vector(Na);

    for (int i = 0; i < Na; ++i) {
        null_vector[i] = svd.v[i][Na-1];
    }

    utils::print(null_vector, "Unit vector in N(A)");
   

    // 3. Use the Singular Value Decomposition to compute the solution x to Ax = b. State the solution x.
    VecDoub x(Na);
    svd.solve(b, x);
    utils::print(x, "solution x");

    // 4. State an estimate of the accuracy on the solution x. State an explanation of how you computed the accuracy.
    VecDoub error_est(Na);
    double sum;
    for (int j = 0; j < Na; j++) {
        sum = 0;
        for (int i = 0; i < Na; i++) {
        if (svd.w[i] > 1e-15) { // to avoid dividing by zero
            sum += pow(svd.v[j][i] / svd.w[i], 2);
        }
        }
        error_est[j] = sqrt(sum);
    }
    utils::print(error_est, "Estimate of accuracy on x");

    return 0;

}
