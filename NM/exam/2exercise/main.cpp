#include <iostream>
#include "nr3.h"
#include "utils.h"
#include "ludcmp.h"
#include <cmath>
#include "qrdcmp.h"
#include "roots_multidim.h"

Doub norm2(VecDoub x) {
  Doub res = 0;
  for (int i = 0; i < x.size(); i++) {
    res += x[i] * x[i];
  }
  return sqrt(res);
}

VecDoub equations(VecDoub_I &x)
{
    VecDoub f(4);
    Doub x0 = x[0];
    Doub x1 = x[1];
	Doub x2 = x[2];
	Doub x3 = x[3];
   
	f[0] = 3*x0 + x1*sin(x2) - cos(x0) + cos(x1*x1) + 4.2;
	f[1] = 3*x1 + x0*x2*x3 + sin(x1) - 5.1;
	f[2] = -x1*x1 + x2*(x3*x3) + 3*x2 + 5.2;
	f[3] = x0 + 3*x3 + sin(x2*x2 * x3*x3) + cos(x1) - 2.3;

    return f;
}

// Modified lnsrsch from NR
template <class T>
Doub my_lnsrch(VecDoub_I &xold, const Doub fold, VecDoub_I &g, VecDoub_IO &p,
VecDoub_O &x, Doub &f, const Doub stpmax, Bool &check, T &func) {
	const Doub ALF=1.0e-4, TOLX=numeric_limits<Doub>::epsilon();
	Doub a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
	Doub rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
	Int i,n=xold.size();
	check=false;
	for (i=0;i<n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=0;i<n;i++)
			p[i] *= stpmax/sum;
	for (i=0;i<n;i++)
		slope += g[i]*p[i];
	if (slope >= 0.0) throw("Roundoff problem in lnsrch.");
	test=0.0;
	for (i=0;i<n;i++) {
		temp=abs(p[i])/MAX(abs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
		f=func(x);
		if (alam < alamin) {
			for (i=0;i<n;i++) x[i]=xold[i];
			check=true;
			return alam;
		} else if (f <= fold+ALF*alam*slope) return alam;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(f-fold-slope));
			else {
				rhs1=f-fold-alam*slope;
				rhs2=f2-fold-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc < 0.0) tmplam=0.5*alam;
					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
					else tmplam=-slope/(b+sqrt(disc));
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = f;
		alam=MAX(tmplam,0.1*alam);
	}
}

// Modified newt from NR
template <class T>
void my_newt(VecDoub_IO &x, Bool &check, T &vecfunc, Int max_iter = 200) {
    cout << setw(3) << "k" << setprecision(7)
    << setw(15) << "x0" 
    << setw(15) << "x1"
    << setw(15) << "x2"
    << setw(15) << "x3"
	<< setw(15) << "Backtracking"
    << setw(15) << "lambda"
    << setw(15) << "error est."
    << endl;
    VecDoub diff_prev(x.size(), 1e-8);
	const Int MAXITS=max_iter;
	const Doub TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
	const Doub TOLX=numeric_limits<Doub>::epsilon();
	Int i,j,its,n=x.size();
	Doub den,f,fold,stpmax,sum,temp,test;
	VecDoub g(n),p(n),xold(n);
	MatDoub fjac(n,n);
	NRfmin<T> fmin(vecfunc);
	NRfdjac<T> fdjac(vecfunc);
	VecDoub &fvec=fmin.fvec;
	f=fmin(x);
	test=0.0;
	for (i=0;i<n;i++)
		if (abs(fvec[i]) > test) test=abs(fvec[i]);
	if (test < 0.01*TOLF) {
		check=false;
		return;
	}
	sum=0.0;
	for (i=0;i<n;i++) sum += SQR(x[i]);
	stpmax=STPMX*MAX(sqrt(sum),Doub(n));
	for (its=0;its<MAXITS;its++) {
		fjac=fdjac(x,fvec);
		for (i=0;i<n;i++) {
			sum=0.0;
			for (j=0;j<n;j++) sum += fjac[j][i]*fvec[j];
			g[i]=sum;
		}
		for (i=0;i<n;i++) xold[i]=x[i];
		fold=f;
		for (i=0;i<n;i++) p[i] = -fvec[i];
		LUdcmp alu(fjac);
		alu.solve(p,p);

		Doub lambda = 1.0;
		lambda = my_lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin); // uses modfified lnsrch to track lambda
        
		VecDoub diff = x - xold;
		Doub diff_norm = norm2(diff);
        Doub C = NAN, e = NAN;
        if (its > 0) {
            Doub diff_prev_norm = norm2(diff_prev);
            C = diff_norm / (diff_prev_norm * diff_prev_norm);
            e = C * (diff_norm * diff_norm);
        }

        cout << setw(3) << its + 1 
        << setw(15) << x[0]
        << setw(15) << x[1] 
        << setw(15) << x[2] 
        << setw(15) << x[3] 
		<< setw(15) << ((lambda == 1) ? "No" : "Yes")
		<< setw(15) << lambda
        << setw(15) << e
        << endl;

        diff_prev = diff;

		test=0.0;
		for (i=0;i<n;i++)
			if (abs(fvec[i]) > test) test=abs(fvec[i]);
		if (test < TOLF) {
			check=false;
			return;
		}
		if (check) {
			test=0.0;
			den=MAX(f,0.5*n);
			for (i=0;i<n;i++) {
				temp=abs(g[i])*MAX(abs(x[i]),1.0)/den;
				if (temp > test) test=temp;
			}
			check=(test < TOLMIN);
			return;
		}
		test=0.0;
		for (i=0;i<n;i++) {
			temp=(abs(x[i]-xold[i]))/MAX(abs(x[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX)
			return;
	}
	return;
	throw("MAXITS exceeded in newt");
}

int main(int argc, char const *argv[])
{
	// Exercise 1: With x0 = −0.7, x1 = 1.2, x2 = 2.3, x3 = −4.1 state (with at least 7 digits) the val-
	// ues of the left hand side of the four equations. (HINT: you should get something around
	// (2.36, 6.03, 49.32, −14.12).

    VecDoub x(4);
    x[0] = -0.7;
    x[1] = 1.2;
	x[2] = 2.3;
	x[3] = -4.1;

    VecDoub res = equations(x);

	cout << "F0(x) = " << setprecision(10) << res[0] << endl;
	cout << "F1(x) = " << setprecision(10) << res[1] << endl;
	cout << "F2(x) = " << setprecision(10) << res[2] << endl;
	cout << "F3(x) = " << setprecision(10) << res[3] << endl;

	cout << endl;
	// Exercise 2:
	// Search for a solution to the equations by performing 7 iterations with the globally convergent
	// Newton method with initial guess (x0, x1, x2, x3) = (0, 0, 0, 0). State for each iteration the
	// values of x0, x1, x2, x3 and for each iteration whether backtracking was applied, and if so,
	// what the value of λ was.

    x[0] = 0.0;
    x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;

	bool check;
	my_newt(x, check, equations, 7);

	cout << endl;
	// Exercise 3:
	// State an estimate the accuracy after 7 iterations. The estimate must be based on the data
	// obtained from the 7 iterations and must be stated with a clear argument of how you computed
	// it. Without such an argument, there is no points for the answer
	
	// Estimate taken from print above.

	// If you are reading this, have a nice day :)
    return 0;
}
