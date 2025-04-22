#include <iostream>
#include "nr3.h"
#include "utilities.h"
#include "ludcmp.h"
#include <cmath>
#include "qrdcmp.h"
#include "roots_multidim.h"

const int v = 120;       // kg
const Doub k = 2.5;      // m
const Doub w = 4.0;      // kg/m
const Doub alpha = 2e-7; // kg⁻1

VecDoub equations(VecDoub_I &q)
{
    Doub n = 0.1;
    Doub d = 30.0;

    VecDoub f(8);
    Doub L0 = q[0];
    Doub L = q[1];
    Doub p = q[2];
    Doub x = q[3];
    Doub theta = q[4];
    Doub phi = q[5];
    Doub a = q[6];
    Doub H = q[7];


    f[0] = (a * (cosh(x / a) - 1)) - p;                // p
    f[1] = (2 * a * sinh(x / a)) - L;                  // L
    f[2] = (2 * x + 2 * k * cos(theta)) - d;           // d
    f[3] = (p + k * sin(theta)) - n;                   // n
    f[4] = sinh(x / a) - tan(phi);                     // phi
    f[5] = (1 + v / (w * L0)) * tan(phi) - tan(theta); // theta
    f[6] = L0 * (1 + alpha * H) - L;                   // L
    f[7] = (w * L0) / (2 * sin(phi)) - H;              // H

    return f;
}

int main(int argc, char const *argv[])
{
    VecDoub q(8);
    q[0] = 27;       // L0
    q[1] = 25;       // L
    q[2] = 3;        // p
    q[3] = 14;       // x
    q[4] = M_PI / 8; // theta
    q[5] = M_PI / 16; // phi
    q[6] = 40;       // a
    q[7] = 100;      // H

    bool check;

    newt(q, check, equations);

    printf("L0: %g\n", q[0]);
    printf("L: %g\n", q[1]);
    printf("p: %g\n", q[2]);
    printf("x: %g\n", q[3]);
    printf("theta: %g\n", q[4]);
    printf("phi: %g\n", q[5]);
    printf("a: %g\n", q[6]);
    printf("H: %g\n", q[7]);

    printf("\n");
    printf("L0: %g\n", q[0]);
    printf("H: %g\n", q[7]);


    return 0;
}
