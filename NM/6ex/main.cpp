#include <iostream>
#include "nr3.h"
#include <math.h>
#include <vector>

using namespace std;

double f(double x)
{
    return x - cos(x);
}

void bisection(double x_i, double y_i, double threshold){
    int k = 0;
    double x_ip1, y_ip1, x_im1;
    double dk, dk_m1;

    // Print table header
    std::cout << std::setw(15)  << "k" << " |"
              << std::setw(15) << "x_k" << " |"
              << std::setw(15) << "x_k-x_{k-1}"<< " |"
              << std::setw(15) << "c" << " |"
              << std::setw(15) << "est. of |ε|" << " |"
              << std::endl;
    std::cout << "Results:" << std::setw(75) << "|d_k|  " << std::endl;

    while(true){
        x_ip1 = (x_i + y_i) / 2;

        if( (f(x_ip1) * f(y_i)) < 0){
            y_ip1 = y_i;
        }else{
            y_ip1 = x_i;
        }

        k++;

        // Ensure x_im1 is properly defined before computing c
        dk = (x_ip1 - x_i);
        if (k > 1) {
            std::cout << std::setw(15) << k << " |"
            << std::setw(15) << x_ip1 << " |"
            << std::setw(15) << dk << " |"
            << std::setw(15) << (abs(dk) / abs(dk_m1)) << " |"
            << std::setw(15) << abs(dk) << " |"
            << std::endl;
        } else {
            std::cout << std::setw(15) << k << " |"
            << std::setw(15) << x_ip1 << " |"
            << std::setw(15) << dk << " |"
            << std::setw(15) << "N/A" << " |" // No c value for first iteration
            << std::setw(15) << abs(dk) << " |"
            << std::endl;
        }

        dk_m1 = dk;
        x_im1 = x_i; // Set x_im1 for folowing iteration
        x_i = x_ip1; 
        y_i = y_ip1;

        if (abs(dk) < threshold){
            break;
        }
    }
    std::cout << "Bisection: Solution is " << x_i << std::endl;
}

void regulaFalsi(double x0, double y0, double epsilon)
{
    double xi = x0;
    double yi = y0;
    double fxi, fyi, x_i1, y_i1, dk, dk_im1, C, ek;
    int k = 0;

    // Print table header
    std::cout << std::setw(15) << "k" << " |"
              << std::setw(15) << "x_k" << " |"
              << std::setw(15) << "dk: x_k-x_{k-1}" << " |"
              << std::setw(15) << "c" << " |"
              << std::setw(15) << "est. of |ε|" << " |"
              << std::endl;
    std::cout << "Results:" << std::setw(75) << "(-C*dk)/(1-C)" << std::endl;
    do
    {
        fxi = f(xi);
        fyi = f(yi);
        x_i1 = xi - ((xi - yi) / (fxi - fyi)) * fxi;

        if (f(x_i1) * fyi < 0)
        {
            y_i1 = yi;
        }
        else
        {
            y_i1 = xi;
        }
        k++;

        dk_im1 = dk;
        dk = x_i1 - xi;
        yi = y_i1;
        xi = x_i1;
        C = dk / dk_im1;
        ek = (-C / (1 - C)) * dk;
        if (k > 1)
        {
            std::cout << std::setw(15) << k << " |"
                      << std::setw(15) << x_i1 << " |"
                      << std::setw(15) << dk << " |"
                      << std::setw(15) << C << " |"
                      << std::setw(15) << ek << " |"
                      << std::endl;
        }
        else
        {
            std::cout << std::setw(15) << k << " |"
                      << std::setw(15) << x_i1 << " |"
                      << std::setw(15) << dk << " |"
                      << std::setw(15) << "N/A" << " |" // No c value for the first iteration
                      << std::setw(15) << "N/A" << " |" // No c value for the first iteration
                      << std::endl;
        }

    } while (abs(dk) > epsilon);
}

int sign(double x){
    if (x<0){
        return -1;
    }
    return 1;
}

void ridders(double x0, double y0, double epsilon){
    double xi = x0;
    double yi = y0;
    double zi, x_i1, y_i1, dk, dk_m1, C, ek;
    int k = 0;

    cout << "RIDDERS" << endl;
    do
    {
        zi = (xi + yi) / 2;
        x_i1 = zi + (zi - xi) * (sign(f(xi)-f(yi))*f(zi) / sqrt(f(zi) * f(zi) - f(xi) * f(yi)));
        if (f(x_i1) * f(zi) < 0) {
            y_i1 = zi;
        } else if (f(x_i1) * f(yi) < 0) {
            y_i1 = yi;
        } else {
            y_i1 = xi;
        }
        dk_m1 = dk;
        dk = x_i1 - xi;
        C = abs(dk) / pow(abs(dk_m1),3);
        xi = x_i1;
        yi = y_i1;

        ek = C * pow(abs(dk),3);
        std::cout << std::setw(15) << k << " |"
              << std::setw(15) << x_i1 << " |"
              << std::setw(15) << dk << " |"
              << std::setw(15) << C << " |"
              << std::setw(15) << ek << " |"
              << std::endl;
        k++;
        } while (abs(dk) > epsilon);
    
}

int main(int argc, char const *argv[])
{
    bisection(0, M_PI/2, 1e-8);

    regulaFalsi(0, M_PI / 2, 1e-16); // 1e-16 er nede på et niveau hvor machine precision er problematisk, så 1e-8 ville være bedre at bruge i virkelighed ish

    ridders(0, M_PI/2, 1e-8);

    return 0;
}
