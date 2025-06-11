#include <cassert>
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
// #include <print>
#include "nr3.h"
#include "utilities.h"
#include "newt_lekt10.h" 

#include "lektion10_ODE.h"



VecDoub vecFunc(double x, VecDoub y){
// Denne funktion indeholder mat funktionerne 
    VecDoub f(y.size());
    
    f[0] = y[0]*y[1];       // du/dx
    f[1] = -(y[0]*y[0]);    // dv/dx

    return f;
}



/*
Pga. jeres problemer med at implementere Trapez-metoden i Mandatory 4 får I her et testproblem som er hugget fra kemisk kinetik:

y1'(t)= - 0.05 y1(t) + 10000 y2(t)y3(t); y1(0)=1
y2'(t)= 0.05 y1(t) - 10000 y2(t)y3(t) - 30000000y2(t)^2; y2(0)=0
y3'(t)=30000000y2(t)^2; y3(0)=0

Brug Trapez-metoden med h=0.005 og find y1(5),y2(5),y3(5) (altså svarende til 1000 skridt). 
I bør få ca. {0.8713, 0.000022, 0.1286} (medmindre jeg har lavet en fejl hvilket man aldrig kan udelukke :)). 
Hvis I ikke har implementeret Trapezmetoden korrekt bliver dette meget tydeligt når I forsøger at løse dette problem.

*/

VecDoub testAfTrapz(double x, VecDoub y){
    VecDoub f(y.size());
 
    // y1 == y[0]

    f[0] = -0.05 * y[0] + 10000 * y[1] * y[2]; // y1'
    f[1] = 0.05 * y[0] - 10000 * y[1] * y[2] - 30000000 * y[1]*y[1]; // y2'
    f[2] = 30000000 * pow(y[1],2); // y3'

    return f;
}

VecDoub handin4(double t, VecDoub v){
// Denne funktion indeholder mat funktionerne 
    VecDoub f(v.size());
    
    f[0] = exp(-t)*cos(v[1])+pow(v[2],2)-v[0];    // du/dx
    f[1] = cos(pow(v[2],2))-v[1];  // dv/dx
    f[2] = cos(t)*exp(-pow(v[0],2))-v[2];

    return f;
}



int main(int, char**){  
/*
u'(x) = u(x) v(x)       u(0) = 1 
v'(x) = -u^2(x)         v(0) = 1
0 <= x  <= 10       for n = 5, 10, 20, 40
*/
    VecDoub y(2);
    y[0] = 1;   // initial condition
    y[1] = 1;   // -------||--------

    // explicit metoder
    richardsonFunc(eulerODE, vecFunc, y, 0, 10, 1, "Euler ODE"); // order er 1 order

    // richardsonFunc(midpointODE, vecFunc, y, 0, 10, 2, "Midpoint ODE"); // er 2 order --> slide 4 lektion 10

    
    // richardsonFunc(RK4, vecFunc, y, 0, 10, 4, "Runge-kutta 4"); // er 4 order --> slide 4 lektion 10

    // implicit metoder.
    // richardsonFunc(trapazODE, vecFunc, y, 0, 10, 2, "trapz");


    // VecDoub ytest(3);
    // ytest[0] = 1; 
    // ytest[1] = 0;
    // ytest[2] = 0;
    // richardsonFunc(trapazODE, testAfTrapz, ytest, 0, 5, 2, "trapz");

    // VecDoub yhandin4(3);
    // yhandin4[0] = 1;   // initial condition
    // yhandin4[1] = 2;   // -------||--------
    // yhandin4[2] = 3;   // -------||--------
    // richardsonFunc(trapazODE, handin4, yhandin4, 0, 5, "trapz");
    
    
}
