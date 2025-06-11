#include <iostream>
#include <cmath>
#include "nr3.h"
#include "quadrature.h"
#include "derule.h"

using namespace std;

double f1(double x, double del)
{
    return cos(x * x) * exp(-x);
}

double f2(double x, double delta)
{
    return sqrt(x) * cos(x * x) * exp(-x);
}

double f3(double x, double del)
{
    return (1 / sqrt(x)) * cos(x * x) * exp(-x);
}

double f4(double x, double delta)
{
    return 1000 * exp(-1.0 / x) * exp(-1.0 / (1.0 - x));
}

void derule_table(double (*f)(double, double), double a, double b, double hmaxx = 4.3, double tol = 1e-16)
{
    DErule<Doub (*)(Doub, Doub)> de(f, a, b, hmaxx);

    vector<double> A;
    int i = 0;
    double curr = 0, diff = 0;

    cout << setw(3) << "i"
         << setw(15) << "A(h_i)"
         << setw(25) << "A(h_(i-1)) - A(h_i)" << endl;

    while (true)
    {
        curr = de.next();
        A.push_back(curr);

        cout << setw(3) << i + 1
             << setw(15) << curr;

        if (i == 0)
            cout << setw(25) << "-";
        else
        {
            diff = A[i - 1] - A[i];
            cout << setw(25) << diff;

            if (abs(diff) < tol)
                break;
        }

        cout << endl;
        i++;
    }

    cout << "\nFinal estimate: " << A.back() << endl;
}



int main(int argc, char const *argv[])
{
    cout << "\n------- f1 -------" << endl;
    derule_table(f1, 0, 1, 4.3);

    cout << "\n------- f2 -------" << endl;
    derule_table(f2, 0, 1, 4.3);

    cout << "\n------- f3 -------" << endl;
    derule_table(f3, 0, 1, 4.3, 1e-15); // Ved tol = 1e-16 stopper den ikke

    cout << "\n------- f4 -------" << endl;
    derule_table(f4, 0, 1, 4.3, 1e-14);

    return 0;
}