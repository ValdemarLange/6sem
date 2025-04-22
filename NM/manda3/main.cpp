#include <iostream>
#include <cmath>
#include "nr3.h"
#include "quadrature.h"


double f(double x)
{
    return (cos(x * x * x) * exp(-x)) / sqrt(x);
}

double fde(double x, double del)
{
    return (cos(x * x * x) * exp(-x)) / sqrt(x);
}

double midPoint(double (*f)(double), double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += f(a + h * (i + 0.5));
    }
    return h * sum;
}

// Test function for different values of N
void testN(double (*f)(double), double a, double b, double e)
{
    double xold = midPoint(f, a, b, 2);
    int f_calc = 2;
    double richerror = 1e9;
    int i = 0;
    double x = 0;
    int N = 0;
    double diffold = 0;
    double diff = 0;
    double h = (b - a) / 2;

    std::cout << "     N           h            A(h_i)     A(h_(i-1))-A(h_i)  Rich. error     alp^k     f_computation\n";
    std::cout << "------------------------------------------------------------------------------------\n";
    
    while (richerror > e)
    {
        i++;
        N = std::pow(2, i - 1) + 1; // 2, 3, 5, 9,..
        h = (b - a) / N;
        x = midPoint(f, a, b, N);
        f_calc += N;
        diffold = diff;
        diff = xold - x;
        xold = x;

        if (i > 2)
        {
            richerror = std::abs(diff) / (std::pow(2.0, 2) - 1.0);
        }

        std::cout << std::setw(8) << N
                  << std::setw(15) << h 
                  << std::setw(15) << x;
                  
        if (i > 1){
            std::cout << std::setw(15) << diff;
        }
        else{
            std::cout << "        -      ";
        }

        if (i > 2)
        {
            std::cout << std::setw(15) << richerror
                      << std::setw(15) << diff / diffold;
        }
        else
        {
            std::cout << "        -       "
                      << "        -    ";
        }

        std::cout << std::setw(6) << f_calc << "\n";
    }

    std::cout << "\nresult: " << x << "\n";
}



template<class T>
struct DErule : Quadrature {
	Doub a,b,hmax,s;
	T &func;
    int f_counter = 0; // added for counting the number of function evaluations

	DErule(T &funcc, const Doub aa, const Doub bb, const Doub hmaxx=3.7)
		: func(funcc), a(aa), b(bb), hmax(hmaxx) {n=0;}

	Doub next() {
		Doub del,fact,q,sum,t,twoh;
		Int it,j;
		n++;
		if (n == 1) {
			fact=0.25;
            f_counter++; 
			return s=hmax*2.0*(b-a)*fact*func(0.5*(b+a),0.5*(b-a));
		} else {
			for (it=1,j=1;j<n-1;j++) it <<= 1;
			twoh=hmax/it;
			t=0.5*twoh;
			for (sum=0.0,j=0;j<it;j++) {
				q=exp(-2.0*sinh(t));
				del=(b-a)*q/(1.0+q);
				fact=q/SQR(1.0+q)*cosh(t);
				sum += fact*(func(a+del,del)+func(b-del,del));
                f_counter += 2;
				t += twoh;
			}
			return s=0.5*s+(b-a)*twoh*sum;
		}
	}
};


int main(int argc, char const *argv[])
{
    double a = 0;
    double b = 4;

    testN(f, a, b, 1e-3);

    std::cout << "--------------" << std::endl;

    DErule<Doub(Doub, Doub)> derule(fde, a, b, 3.7);
    double curr = derule.next();
    double prev = 0;
    double eps = 1;
    int i = 1;
    printf(" i    A(h_i)    f_computations\n");
    while (abs(eps) > 1e-8)
    {
        printf(" %d: %10g   %d\n", i, curr, derule.f_counter);

        prev = curr;
        curr = derule.next();
        eps = curr - prev;
        i++;

    }
    printf("\nresult: %7g\n", curr);


    return 0;
}
