#include <iostream>
#include "nr3.h"
#include <math.h>



double f(double x){
    return x-cos(x);
}


void bisection(double x_i, double y_i, double threshold){
    int k = 0;
    double x_ip1, y_ip1;
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
        x_i = x_ip1; 
        y_i = y_ip1;

        if(abs(dk) < threshold){
            break;
        }
    }
    std::cout << "Bisection: Solution is " << x_i << std::endl;
}


int sign(double x){
    if(x<0.0){return -1;}
    else{return 1;}
}

void ridder(double x_i, double y_i, double threshold, double k_convergenceOrder){
    int k = 0;
    double x_ip1, y_ip1;
    double z_i;
    double dk, dk_m1;
    double C;

    // Print table header
    std::cout << std::setw(15)  << "k" << " |"
              << std::setw(15) << "x_k" << " |"
              << std::setw(15) << "x_k-x_{k-1}" << " |"
              << std::setw(15) << "c" << " |"
              << std::setw(15) << "est. of |ε|" << " |"
              << std::endl;
    std::cout << "Results:" << std::setw(75) << "C*|d_k|^3" << std::endl;

    while(true){
        z_i = (x_i + y_i) / 2;
        x_ip1 = z_i + (z_i - x_i) * ( sign(f(x_i)-f(y_i)) * f(z_i) ) / ( sqrt( pow(f(z_i),2) - f(x_i) * f(y_i) ) );

        if((f(x_ip1) * f(z_i)) < 0){
            y_ip1 = z_i;
        }else if((f(x_ip1) * f(y_i)) < 0){
            y_ip1 = y_i;
        }else{
            y_ip1 = x_i;
        }

        k++;

        dk = (x_ip1 - x_i);
        if (k > 1) {
            C = abs(dk) / pow(abs(dk_m1),k_convergenceOrder);
            std::cout << std::setw(15) << k << " |"
            << std::setw(15) << x_ip1 << " |"
            << std::setw(15) << dk << " |"
            << std::setw(15) << C << " |"
            << std::setw(15) << (C * pow(abs(dk),k_convergenceOrder)) << " |"
            << std::endl;
        } else {
            std::cout << std::setw(15) << k << " |"
            << std::setw(15) << x_ip1 << " |"
            << std::setw(15) << dk << " |"
            << std::setw(15) << "N/A" << " |" // No c value for the first iteration
            << std::setw(15) << "N/A" << " |" // No c value for the first iteration
            << std::endl;
        }

        dk_m1 = dk;
        x_i = x_ip1; 
        y_i = y_ip1;

        if(abs(dk) < threshold){
            break;
        }
    }
    std::cout << "Ridder: Solution is " << x_i << std::endl;
}


void regulaFalsi(double x_i, double y_i, double threshold){
    int k = 0;
    double x_ip1, y_ip1;
    double dk, dk_m1;

    // Print table header
    std::cout << std::setw(15)  << "k" << " |"
              << std::setw(15) << "x_k" << " |"
              << std::setw(15) << "x_k-x_{k-1}" << " |"
              << std::setw(15) << "c" << " |"
              << std::setw(20) << "est. of |ε|" << " |"
              << std::endl;
    std::cout << "Results:" << std::setw(80) << "-C/(1-C) * |d_k|" << std::endl;

    while(true){
        x_ip1 = x_i - ((x_i - y_i) / (f(x_i) - f(y_i))) * f(x_i);
        
        if (f(x_ip1)*f(y_i)){
            y_ip1 = y_i;
        }else{
            y_ip1 = x_i;
        }
        
        k++;

        dk = (x_ip1 - x_i);
        if (k > 1) {
            double C = (abs(dk) / abs(dk_m1)); // Ensure dk_m1 is properly defined before computing c
            std::cout << std::setw(15) << k << " |"
            << std::setw(15) << x_ip1 << " |"
            << std::setw(15) << dk << " |"
            << std::setw(15) << C << " |"
            << std::setw(20) << (( (-C) / (1-C)) * dk) << " |"
            << std::endl;
        } else {
            std::cout << std::setw(15) << k << " |"
            << std::setw(15) << x_ip1 << " |"
            << std::setw(15) << dk << " |"
            << std::setw(15) << "N/A" << " |" // No c value for the first iteration
            << std::setw(20) << "N/A" << " |" // No c value for the first iteration
            << std::endl;
        }

        dk_m1 = dk;
        x_i = x_ip1; 
        y_i = y_ip1;

        if(abs(dk) < threshold){
            break;
        }
    }
    std::cout << "Regula Falsi: Solution is " << x_i << std::endl;

}


double derivative_of_f(double (*f_func)(double), double x, double h = 1e-6){
    return (f_func(x + h) - f_func(x)) / h;
}


void newton(double x_i, double threshold, double k_convergenceOrder){
    int k = 0;
    double x_ip1;
    double dk, dk_m1;

    // Print table header
    std::cout << std::setw(15)  << "k" << " |"
    << std::setw(15) << "x_k" << " |"
    << std::setw(15) << "x_k-x_{k-1}" << " |"
    << std::setw(15) << "c" << " |"
    << std::setw(20) << "est. of |ε|" << " |"
    << std::endl;
    std::cout << "Results:" << std::setw(80) << "C * |d_k|^2" << std::endl;

    while(true){
        x_ip1 = x_i - (1 / derivative_of_f(f, x_i)) * f(x_i);
        k++;

        dk = (x_i - x_ip1 );
        if (k > 1) {
            double C = (abs(dk) / pow(abs(dk_m1), k_convergenceOrder) );
            std::cout << std::setw(15) << k << " |"
            << std::setw(15) << x_ip1 << " |"
            << std::setw(15) << dk << " |"
            << std::setw(15) << C << " |"
            << std::setw(20) << (C * pow(abs(dk), k_convergenceOrder)) << " |"
            << std::endl;
        } else {
            std::cout << std::setw(15) << k << " |"
            << std::setw(15) << x_ip1 << " |"
            << std::setw(15) << dk << " |"
            << std::setw(15) << "N/A" << " |" // No c value for the first iteration
            << std::setw(20) << "N/A" << " |" // No c value for the first iteration
            << std::endl;
        }

        dk_m1 = dk;
        x_i = x_ip1; 

        if(abs(dk) < threshold){
            break;
        }
    }
    std::cout << "Newton: Solution is: " << x_i << std::endl;
}


void secant(double x_i, double x_im1, double threshold, double k_convergenceOrder){
    int k = 0;
    double x_ip1;
    double dk, dk_m1;

    // Print table header
    std::cout << std::setw(15)  << "k" << " |"
    << std::setw(15) << "x_k" << " |"
    << std::setw(15) << "x_k-x_{k-1}" << " |"
    << std::setw(15) << "c" << " |"
    << std::setw(20) << "est. of |ε|" << " |"
    << std::endl;
    std::cout << "Results:" << std::setw(80) << "C * |d_k|^1.62" << std::endl;

    while(true){
        x_ip1 = x_i - ( (x_i - x_im1) / (f(x_i) - f(x_im1)) ) * f(x_i);
        k++;

        dk = (x_i - x_ip1);
        if (k > 1) {
            double C = (abs(dk) / pow(abs(dk_m1), k_convergenceOrder) );
            std::cout << std::setw(15) << k << " |"
            << std::setw(15) << x_ip1 << " |"
            << std::setw(15) << dk << " |"
            << std::setw(15) << C << " |"
            << std::setw(20) << (C * pow(abs(dk), k_convergenceOrder)) << " |"
            << std::endl;
        } else {
            std::cout << std::setw(15) << k << " |"
            << std::setw(15) << x_ip1 << " |"
            << std::setw(15) << dk << " |"
            << std::setw(15) << "N/A" << " |" // No c value for the first iteration
            << std::setw(20) << "N/A" << " |" // No c value for the first iteration
            << std::endl;
        }

        dk_m1 = dk;
        x_im1 = x_i;
        x_i = x_ip1; 

        if(abs(dk) < threshold){
            break;
        }
    }
    std::cout << "Secant: Solution is: " << x_i << std::endl;
}


int main(){
    //Task 1
    std::cout << "Bisection: Solve f(x) = x - cos(x) = 0: " << std::endl;
    double x_l = 0, x_h = M_PI/2;
    bisection(x_l, x_h, 1e-8);

    //Task 2
    //Ridder
    std::cout << "\n\nRidder: Solve f(x) = x - cos(x) = 0: " << std::endl;
    ridder(x_l, x_h, 1e-8, 3);

    //Regula Falsi
    std::cout << "\n\nRegula Falsi: Solve f(x) = x - cos(x) = 0: " << std::endl;
    regulaFalsi(x_l, x_h, 1e-16);

    //Newton
    std::cout << "\n\nNewton: Solve f(x) = x - cos(x) = 0: " << std::endl;
    newton(x_l, 1e-8, 2);

    //Secant
    std::cout << "\n\nSecant: Solve f(x) = x - cos(x) = 0: " << std::endl;
    secant(x_l, x_h, 1e-8, 1.62);
    
    return 0;
}