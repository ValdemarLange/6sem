#include <iostream>
#include <cmath>
#include "nr3.h"
#include "utils.h"
#include "tridag.h"

using namespace std;

int f_count = 0;
double f(double x, double t){
    f_count++;
    return sin(M_PI * x) *exp(-t);
}

double g(double x){
    return x*x;
}

double a(double t){
    return 0.0;
}

double b(double t){
    return 1.0 + sin(t);
}


double crank_nicolson(int N, double x_target, double t_target){
    double alpha = 4.0;
    
    double h = 1.0 / N;
    double dt = h;
    double dx = h;
    double r = alpha * dt / (dx * dx);

    int M = t_target / dt;           // antal tidsskridt

    VecDoub u_curr(N+1), u_next(N+1);


    for (int i = 0; i < N+1; i++) {
        double x = i * dx;
        u_curr[i] = g(x);  // initialbetingelse
    }
    VecDoub lower(N-1, -r / 2.0);      // subdiagonal
    VecDoub diag(N-1, 1.0 + r);           // hoveddiagonal
    VecDoub upper(N-1, -r / 2.0);     // superdiagonal
    VecDoub rhs(N-1);                 // højreside
    
    for (int n = 0; n < M; n++) {


        double t_next = n * dt;
        for (int j = 1; j < N; j++)
        {
            rhs[j-1] = 0.5*r*u_curr[j-1] + (1.0 - r)*u_curr[j] + 0.5*r*u_curr[j+1]
            + 0.5 * dt * (f(j * dx, t_next+dt) + f(j * dx, t_next));
        }
        
        VecDoub u_solve(N-1); // kun for de indre punkter (dim skal passe :)
        
        // Løs tridiagonal system: A * u_next = rhs
        tridag(lower, diag, upper, rhs, u_solve);
        

        for (int j = 1; j < N; j++)
            u_next[j] = u_solve[j - 1];

        u_next[0] = a(t_next);   // a(t)
        u_next[N] = b(t_next);   // b(t)
        u_curr = u_next;  // opdater
    }

    return u_curr[N/2];
}

int main(){

    vector<double> u_values;
    double alpha = 2.0;  // fordi N fordobles hver iteration. Ikke alpha i crank nicolson

    double x_target = 0.5;
    double t_target = 10.0;

    cout << "\n" << "umi = u(0.5, 10.0)" << endl;


    cout << setw(5) << "i"
        << setw(10) << "N"
        << setw(20) << "umi"
        << setw(20) << "umi-1 - umi"
        << setw(10) << "k"
        << setw(20) << "rich. error est."
        << setw(20) << "f-calls" 
        << endl;

    double err_est = 1.0;
    int i = 0;

    while (abs(err_est) > 1e-7) { 
        int N = 2 * pow(2, i+1); 
        double umi = crank_nicolson(N, x_target, t_target); 

        u_values.push_back(umi);

        cout << setw(5) << i + 1
            << setw(10) << N
            << setw(20) << setprecision(10) << umi;

        if (i >= 1) {
            double diff = u_values[i - 1] - u_values[i];
            cout << setw(20) << diff;

            if (i >= 2) {
                double k = log(abs((u_values[i - 1] - u_values[i - 2]) / diff)) / log(alpha);
                err_est = diff / (pow(alpha, k) - 1);
                cout << setw(10) << setprecision(6) << k
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

        cout << setw(20) << f_count << endl;
        i++;
    }
    cout << endl;

}