#include <iostream>
#include <cassert>
#include <cmath>
#include "nr3.h"
// #include "utilities.h"
#include "utils.h"
#include "tridag.h"

using namespace std;

int f_count = 0;
double f(double x, double t){
    f_count++;
    return x*(1-x) * cos(t) * exp(-t/10);
}

double g(double x){
    return x*x*x*x;
}

double a(double t){
    return 0.0;
}

double b(double t){
    return 1.0;
}

double u(double x, double t){
    if(t == 0){
        return g(x);
    }
    if(x == 0){
        return 0.0; //a(t)
    }
    if(x == 1){
        return 1.0; //b(t)
    }
    throw("u(x,t) not defined");
}


double crank_nicolson(int N, double x_target, double t_target){
    double alpha = 1.0;
    
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


double PDE_parabolic(int N, double alpha, double t_target, 
    double f(double x, double t), double u_func(double x, double t)
){
    // Boundary values
    auto g = [&](double x) { return u_func(x, 0.0); };
    auto a = [&](double t) { return u_func(0.0, t); };
    auto b = [&](double t) { return u_func(1.0, t); };
    
    //Stepsize
    double h = 1.0/(double)(N);
    double dt = h, dx = h;

    double r = alpha * dt / (dx * dx);
    

    VecDoub u_curr(N+1);
    for (int i = 0; i < N+1; i++){
        u_curr[i] = g(i*dx);
    }

    int steps = t_target / dt;
    for (int i = 0; i < steps; i++)
    {
        // Current time
        double t_n = i * dt;

        //Initialize Jacobian
        VecDoub J_a(N-1); //J(i, i-1)
        VecDoub J_b(N-1); //J(i, i)
        VecDoub J_c(N-1); //J(i, i+1)
        VecDoub J_r(N-1);
        for (int j = 1; j < N; j++){
            J_a[j-1] = -1.0/2.0 * r;
            J_b[j-1] = (1+r);
            J_c[j-1] = -1.0/2.0 * r;

            J_r[j-1] = 0.5*r*u_curr[j-1] + (1-r)*u_curr[j] + 0.5*r*u_curr[j+1] + dt/2.0 * (f(j*dx, t_n+dt) + f(j*dx, t_n));

        }


        //Tridag
        VecDoub u(N-1);
        tridag(J_a, J_b, J_c, J_r, u);

        for (int j = 1; j < N; j++){
            u_curr[j] = u[j-1]; 
        }
        u_curr[0] = a(t_n+dt);
        u_curr[N] = b(t_n+dt);

    }
    // std::cout << u_curr[N/2] << std::endl;
    return u_curr[N/2]; //Da vi skal finde u(x,t) => u_curr(N/2, t_target) = u(1/2, 20)
}


void richardson_PDE(double tol, int N, double alpha, double t_target, 
    double f(double x, double t), double u_func(double x, double t), int order) {


        std::cout << std::setw(5) << "Itr:" << " |" 
        << std::setw(10) << "N" << " |"
        << std::setw(15) << "A(hi) " << " |"
        << std::setw(20) << "A(hi-1) - A(hi)" << " |"
        << std::setw(20) << "Rich. alpha^k" << " |"
        << std::setw(20) << "Rich. error est." << " |"
        << std::endl;
        



    double A_h = 0.0, A_h_minus_1 = 0.0;
    double A_diff = 0.0, A_diff_old = 0.0;
    double alpha_k = 0.0;
    double error_est = tol + tol*10;

    int i = 1; //Iteration

    while(abs(error_est) > tol){        
        A_h = PDE_parabolic(N, alpha, t_target, f, u_func);
        A_diff = A_h_minus_1 - A_h;


        std::cout << std::setw(5) << i << " |" 
        << std::setw(10) << std::defaultfloat << N << " |"
        << std::setw(15) << std::fixed << std::setprecision(6) << A_h << " |";
    
        if (i == 1) {
            std::cout << std::setw(20) << "*" << " |"
                << std::setw(20) << "*" << " |"
                << std::setw(20) << "*" << " |"
            << std::endl;
        } else if (i == 2) {
            std::cout << std::setw(20) << std::scientific << std::setprecision(6) << A_diff << " |"
                << std::setw(20) << "*" << " |"
                << std::setw(20) << "*" << " |"
            << std::endl;
        } else {
            alpha_k = A_diff_old / A_diff;
            error_est = (A_h - A_h_minus_1) / (pow(2, order) - 1.0);

            std::cout << std::setw(20) << std::scientific << std::setprecision(6) << A_diff << " |"
                << std::setw(20) << std::fixed << std::setprecision(5) << alpha_k << " |"
                << std::setw(20) << std::defaultfloat << error_est << " |"
            << std::endl;
        }

        //Opdater gamle værdier
        A_h_minus_1 = A_h;
        A_diff_old = A_diff;

        //Opdater iteration og maks itr
        i++;
        N*=2;
        if(i > 25){
            std::cout << "Maks iteration reached" << std::endl;
            break;
        }
    }
}

int main(){

    vector<double> u_values;
    double alpha = 2.0;  // fordi N fordobles hver iteration. Ikke alpha i crank nicolson

    double x_target = 0.5;
    double t_target = 20.0;

    cout << "\n" << "k er jens' order estimate. Hans alp^k er blot 2^k." << endl;
    cout << "\n" << "umi = u(0.5, 20.0)" << endl;


    cout << setw(5) << "i"
        << setw(10) << "N"
        << setw(20) << "umi"
        << setw(20) << "umi-1 - umi"
        << setw(10) << "k"
        << setw(20) << "rich. error est."
        << setw(20) << "f-calls" // F-kald bliver reset efter hver iteration lige nu :)
        << endl;

    double err_est = 1.0;
    int i = 0;

    while (abs(err_est) > 1e-8) { 
        f_count = 0;
        int N = 2 * pow(2, i); // N = 4, 8, 16, ...
        double umi = crank_nicolson(N, x_target, t_target);  // <== her bruger vi CN

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

        cout << setw(20) << f_count << endl;
        i++;
    }


}