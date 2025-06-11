#include <iostream>
#include <cmath>
#include <vector>
#include "nr3.h"
// #include "utilities.h"
#include "utils.h"
#include "tridag.h"


// F(x,y,y')
double F(double y_prime, double y, double x){
    return cos(y*y) + y_prime - exp(x*x);
}

// Fy = dFdy
double Fy(double y_prime, double y, double x){
    return -2 * y * sin(y*y); 
}

// Fy' = dFdy'
double Fy_prime(double y_prime, double y, double x){
    return 1; 
}

VecDoub finite_difference_method( double a, double alpha, double b, double beta, int N,
    VecDoub y,
    double (*F)(double, double, double),
    double (*dFdy)(double, double, double),
    double (*dFdy_prime)(double, double, double) ){

    //Stepsize
    double h = (b-a)/(double)(N);
    
    // Create x-values
    VecDoub x(N+1);
    for(int i = 0; i < N+1; i++){
        x[i] = a + i * h;
    }
    
    //Initialize Jacobian
    VecDoub Ja(N-1); //J(i, i-1)
    VecDoub Jb(N-1); //J(i, i)
    VecDoub Jc(N-1); //J(i, i+1)
    VecDoub phi(N-1);

    //first equation non zero Jacobian elements
    Ja[0] = 0.0;
    Jb[0] = 2.0 + h*h * dFdy( (y[2]-alpha)/(2.0*h), y[1], x[1] );
    Jc[0] = -1.0 + h/2.0 * dFdy_prime( (y[2]-alpha)/(2.0*h) , y[1] , x[1]);
    phi[0] = -alpha + 2.0 * y[1] - y[2] + h*h * F( (y[2] - alpha)/(2.0*h), y[1], x[1]);
    

    //internal equations non zero Jacobian elements
    for(int i = 2; i < N - 1; i++){
        Ja[i-1] = -1.0 - h/2.0 * dFdy_prime( (y[i+1] - y[i-1] )/(2.0*h) , y[i], x[i] );
        Jb[i-1] = 2.0 + h*h * dFdy( (y[i+1] - y[i-1] )/(2.0*h) , y[i], x[i] );
        Jc[i-1] = -1.0 + h/2.0 * dFdy_prime( (y[i+1] - y[i-1] )/(2.0*h) , y[i], x[i] );
        phi[i-1] = -y[i-1] + 2.0*y[i] - y[i+1] + h*h * F( (y[i+1] - y[i-1])/(2.0*h), y[i], x[i]);
    }

    //last equation non zero Jacobian elements
    Ja[N-2] = -1.0 - h/2.0 * dFdy_prime( (beta - y[N-2])/(2.0*h) , y[N-1], x[N-1] );
    Jb[N-2] = 2.0 + h*h * dFdy( (beta - y[N-2])/(2.0*h) , y[N-1], x[N-1] );
    Jc[N-2] = 0.0;
    phi[N-2] = -y[N-2] + 2.0*y[N-1] - beta + h*h * F( (beta-y[N-2])/(2.0*h), y[N-1], x[N-1]);
    
    //Tridag - løser J * delta_y = -phi(y)
    VecDoub delta_y(N-1);
    tridag(Ja, Jb, Jc, -phi, delta_y);
    
    //update y
    for(int i = 1; i < N; i++){
        y[i] += delta_y[i-1];
    }

    return y;
}

void boundary_value_problem_richardson(double tol, double a, double alpha, double b, double beta, int N, int order) {
    std::cout << std::setw(5) << "Itr" 
              << std::setw(10) << "N"
              << std::setw(15) << "h"
              << std::setw(25) << "A(h_i)"
              << std::setw(20) << "A(h_{i-1}) - A(h_i)"
              << std::setw(20) << "Alpha^k"
              << std::setw(20) << "Rich. error"
              << "\n";

    double A_curr = 0.0, A_prev = 0.0, A_diff = 0.0, A_diff_prev = 0.0;
    double alpha_k = 0.0, rich_error = tol * 10;

    std::vector<VecDoub> guesses;
    std::vector<VecDoub> solutions;

    // Første gæt: lineær interpolation
    VecDoub y_init(N + 1); // y0
    y_init[0] = alpha;
    y_init[N] = beta;
    for (int i = 1; i < N; ++i)
        y_init[i] = alpha + (beta - alpha) * i / double(N);

    guesses.push_back(y_init);
    solutions.push_back(finite_difference_method(a, alpha, b, beta, N, y_init, F, Fy, Fy_prime));
    A_curr = solutions[0][N / 2];  // Estimat i midten

    int iteration = 1;

    while (std::abs(rich_error) > tol && iteration < 15) {
        int N_new = N * 2;
        VecDoub refined_guess(N_new + 1); // y_i
        refined_guess[0] = alpha;
        refined_guess[N_new] = beta;

        // Interpolation af tidligere løsning
        for (int i = 1; i < N_new; ++i) {
            if (i % 2 == 0) { // i == lige  (genbrug punkt)
                refined_guess[i] = solutions.back()[i / 2];
            } else {          // i == ulige (gennemsnit af nabopunkter)
                refined_guess[i] = 0.5 * (
                    solutions.back()[(i - 1) / 2] + solutions.back()[(i + 1) / 2]
                );
            }
        }

        VecDoub y_new = finite_difference_method(a, alpha, b, beta, N_new, refined_guess, F, Fy, Fy_prime);
        
        // N_new / 2 er midten af intervallet. Da intervallet er [a, b] er [0, 2] så estimere det her ved y(1) som ønsket i opgaven.
        double A_new = y_new[N_new / 2];

        A_diff = A_prev - A_new;

        std::cout << std::setw(5) << iteration
                  << std::setw(10) << N_new
                  << std::setw(15) << (b - a) / double(N_new)
                  << std::setw(25) << setprecision(10) << A_new;

        if (iteration == 1) {
            std::cout << std::setw(20) << "*" << std::setw(20) << "*" << std::setw(20) << "*" << "\n";
        } else if (iteration == 2) {
            std::cout << std::setw(20) << A_diff
                      << std::setw(20) << "*" << std::setw(20) << "*" << "\n";
        } else {
            alpha_k = A_diff_prev / A_diff;
            rich_error = (A_new - A_prev) / (std::pow(2, order) - 1.0);

            std::cout << std::setw(20) << A_diff
                      << std::setw(20) << alpha_k
                      << std::setw(20) << rich_error << "\n";
        }

        // gem nye værdier og opdater
        guesses.push_back(refined_guess);
        solutions.push_back(y_new);

        A_prev = A_new;
        A_diff_prev = A_diff;
        A_curr = A_new;
        N = N_new;
        iteration++;
    }

    if (iteration >= 25){
        std::cout << "Maximum iterations reached.\n";
    }
    cout << "\nFinal estimate: " << A_curr << "\n";
}


int main(){
    // Exercise 1
    double a = 0.0;
    double alpha = -1.0;
    double b = 2.0;
    double beta = 3.0;

    // husk at sætte til at tage den rigtige værdi. Her y(1) alts midten sÅ N/2
    boundary_value_problem_richardson(1e-10, a, alpha, b, beta, 2, 2);
    // max its sat til at passe med N=32768

    // Exercise 2
    boundary_value_problem_richardson(1e-6, a, alpha, b, beta, 2, 2);
    // Aflæst fra tabellen til N = 2048
    // Med Richardson formlen


}