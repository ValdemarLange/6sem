#include <iostream>
#include "nr3.h"
#include "ludcmp.h"
#include "utilities.h"

using namespace std;
	
int main() {



	// Exercise 1:
	// Solve A x = b using LU decomposition, and print the result.

	MatDoub A(3,3);
	A[0][0] = 1.0;	A[0][1] = 2.0;	A[0][2] = 3.0;
	A[1][0] = 2.0;	A[1][1] = -4.0;	A[1][2] = 6.0;
	A[2][0] = 3.0;	A[2][1] = -9.0;	A[2][2] = -3.0;

	VecDoub b(3);
	b[0] = 5.0;
	b[1] = 18.0;
	b[2] = 6.0;

	util::print(A);

	LUdcmp lu(A);

	VecDoub x(3);
	lu.solve(b, x);

	// evaluate x
	util::print(lu.lu);

	util::print(x);

	// print x

	// cout << "----------------" << endl;
	MatDoub A1(3,3);
	A1[0][0] = 1.0;	A1[0][1] = 3.0;	A1[0][2] = 5.0;
	A1[1][0] = -2.0;	A1[1][1] = 0.0;	A1[1][2] = -1.0;
	A1[2][0] = 2.0;	A1[2][1] = 3.0;	A1[2][2] = 1.0;

	LUdcmp lu1(A1);
	MatDoub L(3,3);
	MatDoub U(3,3);

	lu1.decompose(L, U);

	cout << "---------------- LOWER --------------------" << endl;

	util::print(L,"L");

	cout << "---------------- UPPER --------------------" << endl;

	util::print(U, "U");

	util::print(lu1.lu,"LU");

	for (int i = 0; i < 3; i++)
	{
		cout << lu.indx[i] << endl;
	}
	

	return 0;
}
