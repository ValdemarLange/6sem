#include <iostream>
#include "nr3.h"
#include "utilities.h"
#include "ludcmp.h"
#include <cmath>
#include "qrdcmp.h"
#include "roots_multidim.h"

VecDoub diff(VecDoub x, VecDoub y) {
  VecDoub res = VecDoub(x.size());
  for (int i = 0; i < x.size(); i++) {
    res[i] = x[i] - y[i];
  }
  return res;
}

Doub norm2(VecDoub x) {
  Doub res = 0;
  for (int i = 0; i < x.size(); i++) {
    res += x[i] * x[i];
  }
  return sqrt(res);
}


VecDoub equations(VecDoub_I &x)
{
    VecDoub f(2);
    Doub x0 = x[0];
    Doub x1 = x[1];
   

    f[0] = x0 + 2*sin(x1-x0) - exp(-sin(x1+x0));
    f[1] = x0 * cos(x1) + sin(x0) - 1;

    return f;
}

// Alt det her kan bare bruges fra roots_multidim.h og sÅ ændre navnet pÅ kopieret newt() funktion

// template <class T>
// void lnsrch(VecDoub_I &xold, const Doub fold, VecDoub_I &g, VecDoub_IO &p,
// VecDoub_O &x, Doub &f, const Doub stpmax, Bool &check, T &func) {
// 	const Doub ALF=1.0e-4, TOLX=numeric_limits<Doub>::epsilon();
// 	Doub a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
// 	Doub rhs1,rhs2,slope=0.0,sum=0.0,temp,test,tmplam;
// 	Int i,n=xold.size();
// 	check=false;
// 	for (i=0;i<n;i++) sum += p[i]*p[i];
// 	sum=sqrt(sum);
// 	if (sum > stpmax)
// 		for (i=0;i<n;i++)
// 			p[i] *= stpmax/sum;
// 	for (i=0;i<n;i++)
// 		slope += g[i]*p[i];
// 	if (slope >= 0.0) throw("Roundoff problem in lnsrch.");
// 	test=0.0;
// 	for (i=0;i<n;i++) {
// 		temp=abs(p[i])/MAX(abs(xold[i]),1.0);
// 		if (temp > test) test=temp;
// 	}
// 	alamin=TOLX/test;
// 	alam=1.0;
// 	for (;;) {
// 		for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
// 		f=func(x);
// 		if (alam < alamin) {
// 			for (i=0;i<n;i++) x[i]=xold[i];
// 			check=true;
// 			return;
// 		} else if (f <= fold+ALF*alam*slope) return;
// 		else {
// 			if (alam == 1.0)
// 				tmplam = -slope/(2.0*(f-fold-slope));
// 			else {
// 				rhs1=f-fold-alam*slope;
// 				rhs2=f2-fold-alam2*slope;
// 				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
// 				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
// 				if (a == 0.0) tmplam = -slope/(2.0*b);
// 				else {
// 					disc=b*b-3.0*a*slope;
// 					if (disc < 0.0) tmplam=0.5*alam;
// 					else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
// 					else tmplam=-slope/(b+sqrt(disc));
// 				}
// 				if (tmplam>0.5*alam)
// 					tmplam=0.5*alam;
// 			}
// 		}
// 		alam2=alam;
// 		f2 = f;
// 		alam=MAX(tmplam,0.1*alam);
// 	}
// }
// template <class T>
// struct NRfdjac {
// 	const Doub EPS;
// 	T &func;
// 	NRfdjac(T &funcc) : EPS(1.0e-8),func(funcc) {}
// 	MatDoub operator() (VecDoub_I &x, VecDoub_I &fvec) {
// 		Int n=x.size();
// 		MatDoub df(n,n);
// 		VecDoub xh=x;
// 		for (Int j=0;j<n;j++) {
// 			Doub temp=xh[j];
// 			Doub h=EPS*abs(temp);
// 			if (h == 0.0) h=EPS;
// 			xh[j]=temp+h;
// 			h=xh[j]-temp;
// 			VecDoub f=func(xh);
// 			xh[j]=temp;
// 			for (Int i=0;i<n;i++)
// 				df[i][j]=(f[i]-fvec[i])/h;
// 		}
// 		return df;
// 	}
// };
// template <class T>
// struct NRfmin {
// 	VecDoub fvec;
// 	T &func;
// 	Int n;
// 	NRfmin(T &funcc) : func(funcc){}
// 	Doub operator() (VecDoub_I &x) {
// 		n=x.size();
// 		Doub sum=0;
// 		fvec=func(x);
// 		for (Int i=0;i<n;i++) sum += SQR(fvec[i]);
// 		return 0.5*sum;
// 	}
// };
template <class T>
void my_newt(VecDoub_IO &x, Bool &check, T &vecfunc, Int max_iter = 200) {
    cout << setw(3) << "k"
    << setw(15) << "x0" 
    << setw(15) << "x1"
    << setw(15) << "dx"
    << setw(15) << "c"
    << setw(15) << "e"
    << endl;
    VecDoub dx_old(x.size(), 1e-8);
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
		lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin);
        
        
        VecDoub dx = diff(x, xold);
        Doub dx_norm = norm2(dx);
        Doub C = NAN, e = NAN;
        if (its > 0) {
            Doub dx_old_norm = norm2(dx_old);
            C = dx_norm / (dx_old_norm * dx_old_norm);
            e = C * (dx_norm * dx_norm);
        }


        cout << setw(3) << its + 1 
        << setw(15) << x[0]
        << setw(15) << x[1] 
        << setw(15) << dx_norm
        << setw(15) << C
        << setw(15) << e
        << endl;

        dx_old = dx;

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
	throw("MAXITS exceeded in newt");
}

int main(int argc, char const *argv[])
{
	// Exercise 1:

    VecDoub x(2);
    x[0] = 1;
    x[1] = 1;

    VecDoub fval = equations(x);

    printf("F0(x) = %g\n", fval[0]);
    printf("F1(x) = %g\n", fval[1]);

    cout << "\n ... " << endl;

	// Exercise 2:
	// Newton's method og globally  convergent newton's method (GCN)
	// newt() bruger GCN (justere delta x for at sikre vi gÅr den rigtige vej)


	// exercise 3 og 4: 
    bool check;

    // Startgæt
    x[0] = 1;
    x[1] = 2;

    my_newt(x, check, equations); // modificeret til at stoppe efter 6 iteration samt printe løbende

    printf("x0: %g\n", x[0]);
    printf("x1: %g\n", x[1]);

	// exercise 4: Formel for error kan findes i lek 7 slide 15. Overvej at bruge |d_k| hvis C ikke er pæn eller noget

    // OBS: Jeg har k=1 efter første opdatering hvor x0, x1, dx, c og e er
    // resultatet efter første iteration
    // Jens har k=1 er starttidspunkt, hvor x0 og x1 er startgæt,
    // men dx, c og e hører til det første newton skridt

    // Derudover absolut ingen ide om hvorfor Jens' e har andre tal.

    return 0;
}
