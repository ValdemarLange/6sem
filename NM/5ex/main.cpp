#include <iostream>
#include "nr3.h"
#include <math.h>
#include <vector>

using namespace std;

double f(double x){
    return x - cos(x);
}

void bisection(double xmin, double xmax, double epsilon = 1e-8) {
    double xmid, dx, prev_dx = xmax - xmin;
    int k = 1;
    
    cout << setw(5) << "k"
         << setw(15) << "xmin"
         << setw(15) << "xmax"
         << setw(15) << "dx"
         << setw(10) << "C"
         << setw(20) << "Est. of |ε|" << endl;

    while ((xmax - xmin) > epsilon) {
        xmid = (xmin + xmax) / 2.0;
        dx = xmax - xmin;

        // Udskriv
        double err = abs(dx);
        double C = dx;
        if (k > 1){
            C = dx / prev_dx;
        }

        cout << setw(5) << k
             << setw(15) << xmin
             << setw(15) << xmax
             << setw(15) << dx
             << setw(10) << C
             << setw(20) << err << endl;

        // Bisection opdatering
        if (f(xmin) * f(xmid) < 0)
            xmax = xmid;
        else
            xmin = xmid;

        prev_dx = dx;
        k++;
    }

    cout << "\nResult: " << (xmin + xmax) / 2.0 << endl;
}

void secant(double x0, double x1, double epsilon = 1e-8) {
    double x_prev = x0;
    double x_curr = x1;
    double x_next;
    double x_exact = 0.7390851332; // kendt løsning (til e)
    
    cout << setw(5) << "k"
         << setw(15) << "x_k"
         << setw(15) << "dk"
         << setw(15) << "Est. of C"
         << setw(15) << "Est. of |ε|" << endl;

    int k = 1;
    while (true)
    {
        double f_prev = f(x_prev);
        double f_curr = f(x_curr);
        double x_next = x_curr - (x_curr - x_prev) / (f_curr - f_prev) * f_curr;

        double dk = x_curr - x_next;
        double dk_prev = x_curr - x_prev;
        double C = abs(dk) / pow(abs(dk_prev),1.62);
        double error = C * pow(abs(dk),1.62);

        // Udskriv værdier
        cout << setw(5) << k
             << setw(15) << x_curr
             << setw(15) << dk
             << setw(15) << C
             << setw(15) << error << endl;

        if (abs(x_next - x_curr) < epsilon) {
            cout << "\nResult: " << x_next << endl;
            return;
        }

        k++;
        x_prev = x_curr;
        x_curr = x_next;
    }
}

void regulaFalsi(double x0, double x1, double epsilon = 1e-8) {
    double x = x0;
    double y = x1;
    double x_next;
    int k = 0;
    double dk_prev = 0;

    cout << setw(5) << "k"
         << setw(15) << "x_k"
         << setw(15) << "y_k"
         << setw(15) << "dx"
         << setw(15) << "est. of C"
         << setw(15) << "est. of |ε|" << endl;

    while (true) {
        double fx = f(x);
        double fy = f(y);

        x_next = x - (x - y) / (fx - fy) * fx;

        double dk = x_next - x;
        double C = dk / dk_prev; 
        double error = -C / (1 - C) * dk;

        cout << setw(5) << k
             << setw(15) << x
             << setw(15) << y
             << setw(15) << dk
             << setw(15) << C
             << setw(15) << error << endl;

        if (f(x_next) * fy < 0){
            y = y;
        }
        else {
            y = x;
        }
        x = x_next;
        dk_prev = dk;
        

        // Stopkriterium
        if (abs(x-y) < epsilon) {
            cout << "\nResult: " << x_next << endl;
            return;
        }

        k++;
    }
}

double sign(double val) { // Til brug i ridders
    return (val > 0) - (val < 0);  // returnerer +1, 0 eller -1
}

void ridders(double x0, double x1, double epsilon = 1e-8) {
    double xi = x0;
    double yi = x1;

    if (f(xi) * f(yi) >= 0) {
        cout << "f(x0) and f(x1) must have opposite signs!" << endl;
        return;
    }

    cout << setw(5) << "k"
         << setw(15) << "x_i"
         << setw(15) << "dk"
         << setw(15) << "est. of C"
         << setw(15) << "est. of |ε|" << endl;
         
    int k = 1;
    double dk_prev = 0;
    double x_prev = 0;

    while (k < 100) {
        double zi = (xi + yi) / 2.0;
        double fxi = f(xi);
        double fyi = f(yi);
        double fzi = f(zi);

        double denom = sqrt(fzi * fzi - fxi * fyi);

        double x_next = zi + (zi - xi) * sign(fxi - fyi) * fzi / denom;
        double dk = x_next - x_prev;
        
        double C = abs(dk) / pow(abs(dk_prev), 3);
        double error = C * pow(abs(dk), 3);

        // Terminate
        if (abs(dk) < epsilon) {
            cout << "\nResult: " << x_next << endl;
            return;
        }

        cout << setw(5) << k
             << setw(15) << x_next
             << setw(15) << dk
             << setw(15) << C   // C og e esitmat er noget rod i ridders
             << setw(15) << error << endl; 


        // Opdater interval baseret på fortegn
        if (f(x_next) * f(zi) < 0) {
            xi = zi;
            yi = x_next;
        } else if (f(x_next) * f(yi) < 0) {
            xi = x_next;
        } else {
            yi = xi;
            xi = x_next;
        }
        x_prev = x_next;
        dk_prev = dk;
        k++;
    }
}

double fdiff(double (*f_func)(double), double x, double h = 1e-6){
    return (f_func(x + h) - f_func(x)) / h;
}

void newton(double x0, double epsilon = 1e-8) {
    int k = 0;
    double x = x0;
    double x_next;
    double dk;
    double dk_prev = 0;

    cout << setw(5) << "k"
        << setw(15) << "x_i"
        << setw(15) << "dk"
        << setw(15) << "est. of C"
        << setw(15) << "est. of |ε|" << endl;
        
    while (k < 100) {
        x_next = x - (1 / fdiff(f, x)) * f(x);
        dk = x - x_next;
        double C = abs(dk) / pow(abs(dk_prev), 2);
        double error = C * pow(abs(dk), 2);

        cout << setw(5) << k
             << setw(15) << x
             << setw(15) << dk
             << setw(15) << C 
             << setw(15) << error << endl; 

        // Terminate
        if (abs(dk) < epsilon) {
            cout << "\nResult: " << x_next << endl;
            return;
        }
        dk_prev = dk;
        x = x_next;
        k++;
    }

}

int main(int argc, char const *argv[])
{
    cout << "\n Bisection " << endl;
    bisection(0, M_PI/2, 1e-8);

    cout << "\n Secant " << endl;
    secant(0, M_PI/2, 1e-16);

    cout << "\n False Postion / Regula Falsi " << endl;
    regulaFalsi(0, M_PI/2, 1e-16);

    cout << "\n Ridder" << endl;
    ridders(0, M_PI/2, 1e-16);

    cout << "\n Newton" << endl;
    newton(0, 1e-16);

    return 0;
}
