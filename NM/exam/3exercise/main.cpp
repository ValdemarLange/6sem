#include <iostream>
#include <cmath>
#include "nr3.h"
#include "ludcmp.h"
#include "utils.h"
#include "qrdcmp.h"
#include "roots_multidim.h"

using namespace std;


VecDoub equations(double t, VecDoub x0){
    double x = x0[0];
    double v = x0[1];
    double XF = 250 + 15*t - 5* sqrt(1+t*t);
    double XF_p = 15 - 5*t / sqrt(1+t*t);
    double Ddes = 50 + max(0.0, v * 1.5 + (v * (v - XF_p)) / (2 * 2.0));

    VecDoub eq(2);
    eq[0] = v;  
    eq[1] = 4.0 * (
        1.0 - pow(v / 25.0, 4.0)
        - pow(Ddes / (XF - x), 2.0)  
    );
    return eq;
}

VecDoub midpointMethod(VecDoub (*vecFunc)(double, VecDoub), VecDoub y0, double a, double b, int n, int& f_evals){ 
    VecDoub y = y0;
    VecDoub y_next(y.size()); 

    double x = a;
    double h = (b-a)/n; 

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

vector<VecDoub> midpointPlotData(
    VecDoub (*vecFunc)(double, VecDoub),
    VecDoub x0,
    double a,
    double b
) {
    double h = 0.25;
    vector<VecDoub> data;
    VecDoub y = x0;
    double t = a;

    data.push_back(x0);  // startpoint

    while (t < b ) {
        VecDoub k1 = h * vecFunc(t, y);
        VecDoub k2 = h * vecFunc(t + 0.5*h, y + 0.5*k1);
        y = y + k2;
        t += h;
        data.push_back(y);
    }

    return data;
}


void print_midpoint_table(
    VecDoub y0,
    double a,
    double b)
{
    double tol = 2e-5;
    int max_iters = 20;
    double order = 2;

    std::vector<VecDoub> all_results;
    std::vector<int> f_counts;

    std::string function_names[] = {"x(t)", "x'(t)"};
    int prev_n = 20;
    double prev_val = 0.0;
    double diff = 1e9;

    VecDoub prev_vals(y0.size(), 0.0);
    int f_evals = 0;
    
    int i = 0;
    double rich_err_est = 1e10;
    while (abs(rich_err_est) > tol && i < max_iters) {
        int n = prev_n;

        VecDoub result = midpointMethod(equations, y0, a, b, n, f_evals);

        all_results.push_back(result);
        f_counts.push_back(f_evals);

        if (i >= 2) {
            double delta = all_results[i - 1][0] - all_results[i][0];
            double delta_prev = all_results[i - 2][0] - all_results[i - 1][0];

            double alpha_k = delta_prev / delta;
            rich_err_est = delta / (std::pow(2.0, order) - 1);
        }

        prev_n *= 2;
        i++;
    }

    // print table for each component (x(t), x'(t))
    for (int comp = 0; comp < y0.size(); comp++) {
        std::cout << "Estimate of " << function_names[comp] << " for t = 10" << "\n";
        std::cout << std::setw(5) << "i"
                  << std::setw(10) << "N"
                  << std::setw(20) << "A(h_i)"
                  << std::setw(20) << "A(h_{i-1})-A(h_i)"
                  << std::setw(20) << "alpha^k"
                  << std::setw(20) << "Rich-error"
                  << std::setw(20) << "f-calls"
                  << "\n";

        for (int i = 0; i < all_results.size(); i++) {
            double curr = all_results[i][comp];

            cout << std::setw(5) << i + 1
                 << std::setw(10) << 20 * pow(2,i)
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
                double rich_err = delta / (std::pow(2.0, order) - 1);

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
    // ii)
    VecDoub x0(2);
    x0[0] = 0;
    x0[1] = 15;
    double t0 = -10;

    VecDoub res = equations(t0, x0);
    cout << "\nx''(t0) = "<< setprecision(10) << res[1] << "\n" << endl;

    // iii)
    double a = -10;
    double b = 10;
    int N = 80;
    int f = 0;
    VecDoub three = midpointMethod(equations, x0, a, b, N, f);
    cout << "x(10) = " << three[0] << "\n" << endl;

    ////////// Print t, x, v and XF-x to copy into .csv file and plot with excel :) ///////////////
    vector<VecDoub> result = midpointPlotData(equations, x0, a, b);

    double t = a;
    double h = 0.25;
    cout << "t,x,v/x',XF-x" << endl;
    for (VecDoub data : result) {
        double x = data[0];      // position
        double v = data[1];      // velocity
        double XF = 250 + 15 * t - 5 * sqrt(1 + t*t);  // distance between cars

        cout << t << "," << x << "," << v << "," << XF - x << endl;

        t += h;
    }

    // v)
    print_midpoint_table(x0, a, b);

}


