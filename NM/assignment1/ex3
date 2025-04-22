#include <iostream>
#include <fstream>
#include "nr3.h"
#include "utilities.h"
#include "cholesky.h"
#include "ludcmp.h"
#include <cmath>
#include "svd.h"

//your extra include headers

using namespace std;

int main() {
VecDoub xFilip(82); VecDoub yFilip(82);
ifstream Filip("FilipData.dat");
for(int i = 0; i < 82; i++) {
	Filip >> yFilip[i];
	Filip >> xFilip[i];
}

VecDoub xPont(40); VecDoub yPont(40);
ifstream Pont("PontiusData.dat");
for(int i = 0; i < 40; i++) {
	Pont >> yPont[i];
	Pont >> xPont[i];
}

// your code

//Pontius opgave
std::cout << "PONTIUS OPGAVE" << std::endl;
MatDoub A(40,3);
VecDoub b(40);


for(int i = 0; i < 40; i++){
	A[i][0] = 1;
	A[i][1] = xPont[i];
	A[i][2] = xPont[i]*xPont[i];

	b[i] = yPont[i]/1;
}

//util::print(A);
//util::print(b);



MatDoub C(3,3);
C = util::Transpose(A)*A;

VecDoub c(3);
c = util::Transpose(A)*b;



Cholesky cho(C);
VecDoub a(3);
cho.solve(c,a);

util::print(a, "Cholesky a for pontius");

LUdcmp lu_p(C);
lu_p.solve(c,a);
util::print(a, "LU a for pontius");






std::cout << "\nFILIP OPGAVE" << std::endl;
MatDoub Af(xFilip.size(),11);
VecDoub bf(xFilip.size());

for(int i = 0; i < xFilip.size(); i++){
	for(int j = 0; j < 11; j++){
		Af[i][j] = pow(xFilip[i],j);
	}

	bf[i] = yFilip[i]/1;
}

//util::print(Af, "A filip");
//util::print(bf, "b filip");



MatDoub Cf(11,11);
Cf = util::Transpose(Af)*Af;

VecDoub cf(11);
cf = util::Transpose(Af)*bf;

LUdcmp lu(Cf);
VecDoub aflu(11);
lu.solve(cf,aflu);
util::print(aflu, "LU filip");
std::cout << std::endl;

	// Ikke positiv definit derfor virker cholesky ikke
// std::cout << "Cholesky" << std::endl;
// Cholesky chof(Cf);
// VecDoub af(11);
// cho.solve(cf,af);

// ----------------------------- SVD --------------------
cout << "---------------------- SVD _--------------------" << endl;
// ------------ SVD Pont -----------------
SVD svdp(A);
VecDoub xp(3);
svdp.solve(b, xp);

util::print(xp, "Pont SVD x:");


/// ---------- SVD Filip -----------------------------
int Mf = 82;
int Nf = 11;

SVD svd(Af);
VecDoub x(Nf);
svd.solve(bf, x, 1e-15);

util::print(x, "Filip SVD x:");
cout << "rank: " << svd.rank() << endl;
cout << "nullity " << svd.nullity() << endl;
util::print(svd.range(), "Filip range: ");
util::print(svd.nullspace(), "Filip nullspace: ");



return 0;
}