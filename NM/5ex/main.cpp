#include <iostream>
#include "nr3.h"
#include <math.h>
#include <vector>

using namespace std;

double f(double x){
    return x - cos(x);
}

vector<vector<double>> bisection(double x0, double y0, double epsilon = 1e-8){
    vector<double> yk = {y0};
    vector<double> xk = {x0};
    vector<double> dk = {};
    int k = 0;
    double interval = abs(yk[k] - xk[k]);
    while (interval > epsilon){
        xk.push_back( (xk[k] + yk[k]) / 2 );

        if( f(xk[k+1])*f(yk[k]) < 0){
            yk.push_back( yk[k] );
        }
        else{
            yk.push_back( xk[k] );
        }
        interval = abs(yk[k] - xk[k]);
        dk.push_back( xk[k] - xk[k-1] );
        k++;
    }
    vector<vector<double>> out = {xk, yk, dk};
    return out;
}


int main(int argc, char const *argv[])
{
    cout << "\n Bisection " << endl;
    
    vector<vector<double>> res =  bisection(0, M_PI/2);
    cout << "k | xk | dk | yk | C" << endl;
    for (int k = 1; k < res[1].size(); k++)
    {
        cout << setw(4) << k << " | " << setw(8) << res[0][k] << " | " << setw(12) << res[2][k] << " | " << setw(8) << res[1][k] << " | " << res[2][k] / res[2][k-1] << endl;
    }
    

    return 0;
}
