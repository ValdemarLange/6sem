#include <iostream>
#include <cmath>
#include "nr3.h"
#include "quadrature.h"
#include "derule.h"

double f1(double x, double del)
{
    return cos(x * x) * exp(-x);
}

double f2(double x, double delta)
{
    return sqrt(x) * exp(-x);
}

double f3(double x, double del)
{
    return 1 / (sqrt(x))*cos(x * x) * exp(-x);
}

double f4(double x, double delta)
{
    return 1000 * exp(-1 / x) * exp(-1 / (1 - x));
}

int main(int argc, char const *argv[])
{
    // Trapzd<decltype(f1t)> trap(f1t, 0, 1);
    // std::cout << trap.next() << std::endl;

    std::cout << "------- f1 -------" << std::endl;

    DErule<Doub (Doub, Doub)> de1(f1, 0, 1);
    // DErule<T> de1(f1, 0, 1);

    double res = 0; /// LAV OM TIL CURR OG PREV SÅ JEG IKKE KALDER NEXT 2 gange
    double eps = 1;
    int i = 0;
    while (abs(eps) > 1e-8)
    {
        eps = res - de1.next();
        res = de1.next();
        printf("Iteration %d: %g\n", i, res);
        i++;
    }
    std::cout << de1.next() << std::endl;

    std::cout << "------- f2 -------" << std::endl;

    DErule<Doub (Doub, Doub)> de2(f2, 0, 1);

    res = 0;
    eps = 1;
    i = 0;
    while (abs(eps) > 1e-8)
    {
        eps = res - de2.next();
        res = de2.next();
        printf("Iteration %d: %g\n", i, res);
        i++;
    }
    std::cout << de2.next() << std::endl;

    std::cout << "------- f3 -------" << std::endl;

    DErule<Doub (Doub, Doub)> de(f3, 0, 1, 3.7);

    res = 0;
    eps = 1;
    i = 0;
    while (abs(eps) > 1e-8)
    {
        eps = res - de.next();
        res = de.next();
        printf("Iteration %d: %g\n", i, res);
        i++;
    }
    std::cout << de.next() << std::endl;

    std::cout << "------- f4 -------" << std::endl;

    DErule<Doub (Doub, Doub)> de4(f4, 0, 1);

    res = 0;
    eps = 1;
    i = 0;
    while (abs(eps) > 1e-8)
    {
        eps = res - de4.next();
        res = de4.next();
        printf("Iteration %d: %g\n", i, res);
        i++;
    }
    std::cout << de4.next() << std::endl;

    return 0;
}