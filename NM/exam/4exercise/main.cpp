#include <iostream>
#include <cmath>
#include "nr3.h"
#include "utils.h"

using namespace std;

double f(double x)
{
    return exp(x*x*x) * sqrt(x*(2-x));
}


// Simpson's rule // n skal være et ulige tal
double simpson(double (*f)(double), double a, double b, int n)
{
    double h = (b - a) / (n-1);
    double sum = (f(a) + f(b));
    for (int i = 1; i < n-1; i++)
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
    vector<double> h_vec; // Gemmer h værdier

    cout << setw(3) << "i"
         << setw(15) << "N-1"
         << setw(15) << "A(h_i)"
         << setw(20) << "A(h_(i-1))-A(h_i)"
         << setw(15) << "alpha^k"
         << setw(15) << "Rich-error"
         << setw(15) << "f-calls"
         << setw(15) << "order k" 
         << setw(15) << "h" << endl;

    for (int i = 1; i < max_iter+1; i++) {
        int N_m1 = pow(2, i);
        int N = N_m1 +1;
        double h = (b-a) / N_m1;

        double A_i = simpson(f, a, b, N);
        A.push_back(A_i);
        h_vec.push_back(h);

        cout << setw(3) << i
             << setw(15) << N_m1
             << setw(15) << A_i;

        if (i >= 2) {
            double diff = A[i - 2] - A[i - 1];
            cout << setw(20) << diff;

            if (i >= 3) {
                double alphak = (A[i - 3] - A[i - 2]) / (A[i - 2] - A[i - 1]);
                double rich = (A[i - 1] - A[i - 2]) / (alphak - 1.0);
                double order = log(abs(alphak)) / log(h_vec[i-2]/h);

                cout << setw(15) << alphak
                     << setw(15) << rich
                     << setw(15) << N+1
                     << setw(15) << order;
            }
        }
        cout << setw(15) << h << endl;
    }

    cout << "\nFinal estimate: " << setprecision(10) << A.back() << endl;

}


int main(int argc, char const *argv[])
{
    // Exercise 1. With N −1 = 2k; k = 1, . . . , 20 use the Simpson Method method to approximate the integral.
    // State the results in a table similar to those used during the course

    double a = 0;
    double b = 2;
    int max_its = 20;

    cout << "\n Approximation of integral using Simpson" << endl;
    simpson_table(f, a, b, max_its);    

    cout << endl;

    return 0;
}
