#include <iostream>
#include <cmath>
#include "nr3.h"
#include "utilities.h"


double f(double x, double a)
{
    return ( cos(x) * exp(-x*x*x) ) / (pow((x - a),3.0/4.9));
}

// Extended Midpoint (Rectangle) method
double midpoint(double (*f)(double, double), double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += f(a + h * (i + 0.5), a);
    }
    return h * sum;
}
void midpoint_table(double (*f)(double, double), double a, double b, int max_iter = 16) {
    vector<double> A;  // Gemmer A(h_i)
    int n = 1;
    n = 3;
    double h = 2;
    double diff = 0;
    double alp_k = 0;
    double rich = 0;
    double order = 0;

    cout << setw(3) << "i"
         << setw(10) << "N"
         << setw(15) << "A(h_i)"
         << setw(20) << "A(h_(i-1)) - A(h_i)"
         << setw(15) << "alp^k"
         << setw(15) << "Rich-error"
         << setw(15) << "f-calls"
         << setw(15) << "order k" << endl;

    for (int i = 0; i < max_iter; i++) {
        double A_i = midpoint(f, a, b, n);
        A.push_back(A_i);

        if (i > 0) {
            diff = A[i - 1] - A[i];
        }
        if (i > 1) {
            alp_k = (A[i - 2] - A[i - 1]) / (A[i - 1] - A[i]);
            rich = A[i] + (A[i] - A[i - 1]) / (alp_k - 1);
            order = log(abs(alp_k)) / log(h);
        }

        cout << setw(3) << i + 1
             << setw(10) << n
             << setw(15) << A_i;
        if (i > 0) {
            cout << setw(20) << diff;
        }
        if (i > 1) {
             cout << setw(15) << alp_k
             << setw(15) << (A[i] - rich)
             << setw(15) << n
             << setw(15) << order;
        }
        cout << endl;

        n = pow(2, i+2)+1;
    }

    cout << "\nFinal estimate: " << A.back() << endl;
}

// Trapezoidal method
double trapezoid(double (*f)(double), double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));
    for (int i = 1; i < n; i++)
    {
        sum += f(a + h * i);

    }
    return h * sum;
}
void trapezoid_table(double (*f)(double), double a, double b, int max_iter = 16) {
    vector<double> A;  // Gemmer A(h_i)
    int n = 1;
    double h = 2;
    double diff = 0;
    double alp_k = 0;
    double rich = 0;
    double order = 0;

    cout << setw(3) << "i"
         << setw(15) << "A(h_i)"
         << setw(20) << "A(h_(i-1)) - A(h_i)"
         << setw(15) << "alp^k"
         << setw(15) << "Rich-error"
         << setw(15) << "f-calls"
         << setw(15) << "order k" << endl;

    for (int i = 0; i < max_iter; i++) {
        double A_i = trapezoid(f, a, b, n);
        A.push_back(A_i);

        if (i > 0) {
            diff = A[i - 1] - A[i];
        }
        if (i > 1) {
            alp_k = (A[i - 2] - A[i - 1]) / (A[i - 1] - A[i]);
            rich = A[i] + (A[i] - A[i - 1]) / (alp_k - 1);
            order = log(abs(alp_k)) / log(h);
        }

        cout << setw(3) << i + 1
             << setw(15) << A_i;
        if (i > 0) {
            cout << setw(20) << diff;
        }
        if (i > 1) {
             cout << setw(15) << alp_k
             << setw(15) << (A[i] - rich)
             << setw(15) << n+1
             << setw(15) << order;
        }
        cout << endl;

        n *= h;
    }

    cout << "\nFinal estimate: " << A.back() << endl;
}
// Simpson's rule // n skal være et lige tal
double simpson(double (*f)(double), double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = (f(a) + f(b));
    for (int i = 1; i < n; i++)
    {
        if (i % 2 == 1){
        sum += 4.0 *f(a + h * i);
        } else {
        sum += 2.0 *f(a + h * i);
        }
    }
    return (h/3.0) * sum;
}
void simpson_table(double (*f)(double), double a, double b, int max_iter = 16) {
    vector<double> A;  // Gemmer A(h_i)
    int n = 1;
    double h = 2;
    double diff = 0;
    double alp_k = 0;
    double rich = 0;
    double order = 0;

    cout << setw(3) << "i"
         << setw(15) << "A(h_i)"
         << setw(20) << "A(h_(i-1)) - A(h_i)"
         << setw(15) << "alp^k"
         << setw(15) << "Rich-error"
         << setw(15) << "f-calls"
         << setw(15) << "order k" << endl;

    for (int i = 0; i < max_iter; i++) {
        double A_i = simpson(f, a, b, n);
        A.push_back(A_i);

        if (i > 0) {
            diff = A[i - 1] - A[i];
        }
        if (i > 1) {
            alp_k = (A[i - 2] - A[i - 1]) / (A[i - 1] - A[i]);
            rich = A[i] + (A[i] - A[i - 1]) / (alp_k - 1);
            order = log(abs(alp_k)) / log(h);
        }

        cout << setw(3) << i + 1
             << setw(15) << A_i;
        if (i > 0) {
            cout << setw(20) << diff;
        }
        if (i > 1) {
             cout << setw(15) << alp_k
             << setw(15) << (A[i] - rich)
             << setw(15) << (n+1)
             << setw(15) << order;
        }
        cout << endl;

        n *= h;
    }

    cout << "\nFinal estimate: " << A.back() << endl;
}


int main(int argc, char const *argv[])
{
    double a = 1;
    double b = 5;

    cout << "\nMidpoint" << endl;
    midpoint_table(f, 0, 1, 20);
    
    // cout << "\nTrapezoidal" << endl;
    // trapezoid_table(f, 0, 1);

    // cout << "\nSimpson" << endl;
    // simpson_table(f, 0, 1, 10);

    return 0;
}
