#include <iostream>
#include <cassert>
#include <cmath>
#include "nr3.h"
// #include "utilities.h"
#include "utils.h"
#include "tridag.h"


double f(double x, double t){
    return x*(1-x) * cos(t) * exp(-t/10);
}

double g(double x){
    return x*x*x*x;
}

double u_func(double x, double t){
    if(t == 0){
        return g(x);
    }
    if(x == 0){
        return 0.0; //a(t)
    }
    if(x == 1){
        return 1.0; //b(t)
    }
    assert(false);
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


    double alpha = 1.0;
    double t_target = 20.0; //Vi skal finde u(1/2, 20)
    int expected_order = 2;
    int N = 2;
    double tol = 1e-8;

    richardson_PDE(tol, N, alpha, t_target, f, u_func, expected_order);

}