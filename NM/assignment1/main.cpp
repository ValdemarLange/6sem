#include <iostream>
#include <fstream>
#include "nr3.h"
#include "utilities.h"
// #include "cholesky.h"
// #include "ludcmp.h"
#include <cmath>
#include "svd.h"

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
    util::print(dx, "delta x");
}

int main(int argc, char const *argv[])
{
    ifstream input_A("Ex1A.dat");
    int M, N;
    input_A >> M;
    input_A >> N;
    cout << M << " " << N << endl;
    MatDoub A(M,N);
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            input_A >> A[i][j];
        }
    }
    

    ifstream input_b("Ex1b.dat");
    int Mb, Nb;
    input_b >> Mb;
    input_b >> Nb;
    cout << Mb << " " << Nb << endl;
    VecDoub b(Mb);
    for (int i = 0; i < Mb; i++)
    {
        input_b >> b[i];
    }    

    SVD svd(A);
    util::print(svd.w, "W");

    VecDoub x(N);
    svd.solve(b,x);
    util::print(x, "solution x");

    print_error(svd);
    
    VecDoub r(N);
    r = A*x;
    for (int i = 0; i < r.size(); i++)
    {
        r[i] = r[i] - b[i];
    }
    util::print(r, "r");

    VecDoub sigma(M);
    double delta = 1;
    for (int i = 0; i < sigma.size(); i++)
    {
        sigma[i] = max(delta, abs(r[i]));
    }
    util::print(sigma, "sigma_i");

    MatDoub A_weighted(M,N);
    VecDoub b_weighted(Mb);
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A_weighted[i][j] = A[i][j] / sigma[i];
        }
        b_weighted[i] = b[i] / sigma[i];
    }
    cout << "A_0,0: " << A_weighted[0][0] << " B_6: " << b_weighted[6] << endl;


    SVD svd_weighted(A_weighted);

    VecDoub x_weighted(N);
    svd_weighted.solve(b_weighted,x_weighted);
    util::print(x_weighted, "solution x");


    return 0;
}
