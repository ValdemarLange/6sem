#include <iostream>
#include <cmath>
#include "nr3.h"
#include "ludcmp.h"
#include "utils.h"
// #include "utilities.h"
#include "quadrature.h"
#include "qrdcmp.h"
#include "roots_multidim.h"
// #include "newt_lekt10.h"


using namespace std;

VecDoub f10(double x, VecDoub y)
{
    VecDoub f(2);
    f[0] = y[0] * y[1];       // u' = u*v
    f[1] = -y[0] * y[0];      // v' = -u^2
    return f;
}

VecDoub testAfTrapz(double x, VecDoub y){
    VecDoub f(y.size());
 
    // y1 == y[0]

    f[0] = -0.05 * y[0] + 10000 * y[1] * y[2]; // y1'
    f[1] = 0.05 * y[0] - 10000 * y[1] * y[2] - 30000000 * y[1]*y[1]; // y2'
    f[2] = 30000000 * pow(y[1],2); // y3'

    return f;
}

VecDoub eulerODE(VecDoub (*vecFunc)(double, VecDoub), VecDoub y0, double a, double b, int n, int& f_evals){ 
    VecDoub y = y0;
    VecDoub y_next(y.size()); 

    double x = a;
    double h = (b-a)/double(n);

    for (double i = a; i < b; i += h){
       y_next = y + h * vecFunc(x, y);
       f_evals++;
       y = y_next;
       x += h;
    }
    return y_next;
}

VecDoub midpointODE(VecDoub (*vecFunc)(double, VecDoub), VecDoub y0, double a, double b, int n, int& f_evals){ 
    VecDoub y = y0;
    VecDoub y_next(y.size()); 

    double x = a;
    double h = (b-a)/double(n);

    for (double i = a; i < b; i += h){
        VecDoub k1 = h * vecFunc(x,y);
        VecDoub k2 = h * vecFunc(x + 0.5*h, y + 0.5*k1); 
        y_next = y + k2;
        f_evals += 2;
        y = y_next;
        x += h;
    }
    return y_next;
}


// tolf = 1e-7 i roots_multidim.h, ved mindre fejler LUdcmp
int trap_f_count = 0;

VecDoub trapezoidalODE(VecDoub (*vecFunc)(double, VecDoub), VecDoub y0, double a, double b, int n, int& f_evals){
    VecDoub y = y0;
    VecDoub y_next(y0.size());

    double x = a;
    double h = (b-a)/(double)(n);

    for (double i = a; i < b; i += h){
        double x_next = x + h;

        auto phi = [&](const VecDoub& y_guess) -> VecDoub {
            trap_f_count += 2;
            return y_guess - y - (h/2.0) * (vecFunc(x, y) + vecFunc(x_next, y_guess));
        };

        VecDoub guess = y + h * vecFunc(x, y);
        trap_f_count += 1;

        bool check;
        newt(guess, check, phi);
        y = guess;

        x += h;
    }

    return y;
}

VecDoub runge_kutta_4ODE(VecDoub (*vecFunc)(double, VecDoub), VecDoub y0, double a, double b, int n, int& f_evals){ 
    VecDoub y = y0;
    VecDoub y_next(y.size()); 

    double x = a;
    double h = (b-a)/double(n);

    for (double i = a; i < b; i += h){
        VecDoub k1 = h * vecFunc(x,y);
        VecDoub k2 = h * vecFunc(x + 0.5*h, y + 0.5*k1); 
        VecDoub k3 = h * vecFunc(x + 0.5*h, y + 0.5*k2);
        VecDoub k4 = h * vecFunc(x + h, y + k3);
        y_next = y + (1.0/6.0)*k1 + (1.0/3.0)*k2 + (1.0/3.0)*k3 + (1.0/6.0)*k4;
        f_evals += 4;
        y = y_next;
        x += h;
    }
    return y_next;
}

VecDoub leap_frogODE(VecDoub (*vecFunc)(double, VecDoub), VecDoub y0, double a, double b, int n, int& f_evals) {
    double h = (b - a) / double(n); 
    double x = a;

    VecDoub y_prev = y0;  // y₀
    VecDoub y_curr(y0.size());  // y₁
    VecDoub y_next(y0.size());

    // Step 1: Brug Euler som "starter" til at finde y₁
    y_curr = y_prev + h * vecFunc(x, y_prev);
    f_evals++; // f kald til Euler starter
    x += h;

    // Step 2: Leap-Frog iteration
    for (int i = 1; i < n; i++) {
        y_next = y_prev + h * vecFunc(x, y_curr);
            
        f_evals++;

        y_prev = y_curr;
        y_curr = y_next;
        x += h;
    }

    return y_curr;
}


void richardsonFunc(
    VecDoub (*solver)(VecDoub(*)(double, VecDoub), VecDoub, double, double, int, int&), // solver f.eks eulerODE, midpointODE, RK4, trapazODE
    VecDoub (*f)(double, VecDoub),
    VecDoub y0,
    double a,
    double b,
    double method_order,
    const std::string& method_name)
{
    const double tol = 1e-8;
    const int max_iters = 20;

    std::vector<VecDoub> all_results;
    std::vector<int> f_counts;

    std::string function_names[] = {"u(x)", "v(x)", "w(x)"};
    int prev_n = 5;
    // double prev_val = 0.0;
    double diff = 1e9;

    VecDoub prev_vals(y0.size(), 0.0);

    int iter = 0;
    while (diff > tol && iter < max_iters) {
        int n = prev_n;
        int f_evals = 0;
        trap_f_count = 0;

        VecDoub result = solver(f, y0, a, b, n, f_evals);

        all_results.push_back(result);
        if (method_name == "Trapezoidal") {
            f_counts.push_back(trap_f_count); // Tæl f kald for trapezoidal
        } else {
            f_counts.push_back(f_evals); // Tæl f kald for alle andre
        }
        int min_iter_for_diff_check = 0;
        if ( method_name == "Runge-Kutta 4") {
            min_iter_for_diff_check = 1;
        }

        if (iter > min_iter_for_diff_check) {
            // diff = std::abs(result[0] - prev_val); // Tjekker kun konvergens af u(x)
            diff = 0.0;
            for (int i = 0; i < result.size(); i++) { // Tjekker alle komponenter
                double diff_current = std::abs(result[i] - prev_vals[i]);
                if (diff_current > diff)
                    diff = diff_current;
            }
        }

        // prev_val = result[0];
        prev_vals = result;
        prev_n *= 2;
        iter++;
    }

    // Udskriv tabel for hver komponent (u, v, w,...)
    for (int comp = 0; comp < y0.size(); comp++) {
        std::cout << "Method: " << method_name << " | Function: " << function_names[comp] << "\n";
        std::cout << std::setw(5) << "i"
                  << std::setw(20) << "A(h_i)"
                  << std::setw(20) << "A(h_{i-1})-A(h_i)"
                  << std::setw(20) << "alpha^k"
                  << std::setw(20) << "Rich-error"
                  << std::setw(20) << "f-calls"
                  << "\n";

        for (int i = 0; i < all_results.size(); i++) {
            double curr = all_results[i][comp];

            cout << std::setw(5) << i + 1
                 << std::setw(20) << curr;


            if (i == 0) {
                std::cout << std::setw(20) << "*"
                        << std::setw(20) << "*"
                        << std::setw(20) << "*";
            } else if (i == 1) {
                double delta = all_results[i - 1][comp] - curr;
                std::cout << std::setw(20) << delta
                        << std::setw(20) << "*"
                        << std::setw(20) << "*";
            } else {
                double delta = all_results[i - 1][comp] - curr;
                double alpha_k = (all_results[i - 2][comp] - all_results[i - 1][comp]) / delta;
                double rich_err = delta / (std::pow(2.0, method_order) - 1);

                std::cout << std::setw(20) << delta
                        << std::setw(20) << alpha_k
                        << std::setw(20) << rich_err;
            }
            std::cout << std::setw(20) << f_counts[i] << "\n";
            
        }
        std::cout << "\n";
    }
}


int main()
{
    VecDoub y0(2);
    y0[0] = 1.0;  // u(0)
    y0[1] = 1.0;  // v(0)

    // richardsonFunc(eulerODE, f10, y0, 0, 10, 1, "Euler ODE"); // order = 1
    // richardsonFunc(midpointODE, f10, y0, 0, 10, 2, "Midpoint ODE"); // order = 2
    
    // richardsonFunc(trapezoidalODE, f10, y0, 0, 10, 2, "Trapezoidal"); // order = 2 tror jeg. tjek hellere

    // richardsonFunc(runge_kutta_4ODE, f10, y0, 0, 10, 4, "Runge-Kutta 4"); // order = 4

    // richardsonFunc(leap_frogODE, f10, y0, 0, 10, 2, "Leap-Frog"); // order = 2


    /*
    Pga. jeres problemer med at implementere Trapez-metoden i Mandatory 4 får I her et testproblem som er hugget fra kemisk kinetik:

    y1'(t)= - 0.05 y1(t) + 10000 y2(t)y3(t); y1(0)=1
    y2'(t)= 0.05 y1(t) - 10000 y2(t)y3(t) - 30000000y2(t)^2; y2(0)=0
    y3'(t)=30000000y2(t)^2; y3(0)=0

    Brug Trapez-metoden med h=0.005 og find y1(5),y2(5),y3(5) (altså svarende til 1000 skridt). 
    I bør få ca. {0.8713, 0.000022, 0.1286} (medmindre jeg har lavet en fejl hvilket man aldrig kan udelukke :)). 
    Hvis I ikke har implementeret Trapezmetoden korrekt bliver dette meget tydeligt når I forsøger at løse dette problem.

    */

    VecDoub ytest(3);
    ytest[0] = 1; 
    ytest[1] = 0;
    ytest[2] = 0;
    int abe = 0;
    trap_f_count = 0;
    VecDoub res = trapezoidalODE(testAfTrapz, ytest, 0, 5, 1000, abe);
    utils::print(res, "res");
    cout << trap_f_count << endl;
}


