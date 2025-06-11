#include <iostream>
#include <cmath>
#include "nr3.h"
#include "utils.h"

VecDoub vecfunc(VecDoub y)
{
    VecDoub f(2);
    f[0] = y[0] * y[1];
    f[1] = -(y[1] * y[1]);
    return f;
}

VecDoub euler_ODE(VecDoub (*vecfunc)(VecDoub), VecDoub y0, double a, double b, int n)
{
    VecDoub y_current(2);
    y_current[0] = y0[0];
    y_current[1] = y0[1];
    VecDoub y_next(2);

    double h = (b - a) / n;

    for (double i = a; i < b; i += h)
    {
        y_next = y_current + h * vecfunc(y_current);
        y_current = y_next;
    }
    return y_next;
}

int main(int argc, char const *argv[])
{
    VecDoub y0(2);
    y0[0] = 1.0;
    y0[1] = 1.0;
    VecDoub result = euler_ODE(vecfunc, y0, 0, 10, 100);

    std::cout << "u(10) = " << result[0] << ", v(10) = " << result[1] << std::endl;
    return 0;
}
