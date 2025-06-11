#include <cmath>
#include <iostream>

#include "nr3.h"
#include "banded.h"
#include "utils.h"

double f(double x, double y) {
    return 1.0 + x + y;
}

double a0(double x) {
    return 0.0;
}

double a1(double x) {
    return 0.0;
}

double b0(double y) {
    return 0.0;
}

double b1(double y) {
    return 0.0;
}

int inx(int j, int k, int N) {
    return (N - 1) * (k - 1) + (j - 1);
}

double solve_poisson(int N) {
    int n = (N - 1) * (N - 1);
    double h = 1.0 / N;
    double lambda = 0.0;

    // MatDoub A(n, n, 0.0);  // A er n×n, ikke N×N
    VecDoub phi(n, 0.0);

    int m1 = N - 1; // båndbredde ned (underdiagonaler)
    int m2 = N - 1; // båndbredde op (overdiagonaler)
    MatDoub Aband(n, m1 + m2 + 1, 0.0);  // båndmatrix


    for (int j = 1; j < N; ++j) {
    for (int k = 1; k < N; ++k) {
        int i = inx(j, k, N);
        double x = j * h;
        double y = k * h;

        phi[i] = h * h * f(x, y);  // RHS

        // hoveddiagonal
        Aband[i][m1] = 4.0 + h * h * lambda;

        // venstre nabo
        if (j > 1)
            Aband[i][m1 - 1] = -1.0;

        // højre nabo
        if (j < N - 1)
            Aband[i][m1 + 1] = -1.0;

        // nedre nabo
        if (k > 1)
            Aband[i][m1 - (N - 1)] = -1.0;

        // øvre nabo
        if (k < N - 1)
            Aband[i][m1 + (N - 1)] = -1.0;

        // bidrag fra randbetingelser
        if (j == 1)
            phi[i] += a0(y);
        if (j == N - 1)
            phi[i] += a1(y);
        if (k == 1)
            phi[i] += b0(x);
        if (k == N - 1)
            phi[i] += b1(x);
    }
}


    // utils::print(Ab, "Ab");
    // utils::print(phi, "phi");

    Bandec bandecA(Aband, m1, m2);

    VecDoub u(n, 0.0);
    bandecA.solve(phi, u);
    // utils::print(u, "u");

    // cout << "u(0.5, 0.5) = " << u[inx(N / 2, N / 2, N)] << "\n";
    return u[inx(N / 2, N / 2, N)];
}

using namespace std;

int main() {
    vector<double> u_values;
    double alpha = 2.0; // Fordi N fordobles sÅ halveres skridtlængden h hver gang => alpha = h_i-1 / h_i = 2

    cout << "\n" << "umi = u(0.5, 0.5)" << endl;

    cout << setw(5) << "i"
        << setw(10) << "N"
        << setw(20) << "umi"
        << setw(20) << "umi-1 - umi"
        << setw(10) << "k"
        << setw(20) << "rich. error est." << endl;

    double err_est = 1.0;
    int i = 0;
    while (abs(err_est) > 1e-5) { 
        int N = 4 * pow(2, i); // N = 4, 8, 16, ...
        double umi = solve_poisson(N);

        u_values.push_back(umi);

        cout << setw(5) << i + 1
            << setw(10) << N
            << setw(20) << umi;

        if (i >= 1) {
            double diff = u_values[i - 1] - u_values[i];
            cout << setw(20) << diff;

            if (i >= 2) {
                double k = log(abs((u_values[i - 1] - u_values[i - 2]) / diff)) / log(alpha);
                err_est = diff / (pow(alpha, k) - 1);
                cout << setw(10) << k
                    << setw(20) << err_est;
            } else {
                cout << setw(10) << "*"
                    << setw(20) << "*";
            }
        } else {
            cout << setw(20) << "*"
                << setw(10) << "*"
                << setw(20) << "*";
        }
        cout << endl;
        i++;
    }


    return 0;
}
