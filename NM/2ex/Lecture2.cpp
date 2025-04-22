#include <iostream>
#include <fstream>
#include "nr3.h"
#include "utilities.h"
#include "cholesky.h"
#include "ludcmp.h"
// your extra include headers

using namespace std;

int main()
{
	VecDoub xFilip(82);
	VecDoub yFilip(82);
	ifstream Filip("FilipData.dat");
	for (int i = 0; i < 82; i++)
	{
		Filip >> yFilip[i];
		Filip >> xFilip[i];
	}

	VecDoub xPont(40);
	VecDoub yPont(40);
	ifstream Pont("PontiusData.dat");
	for (int i = 0; i < 40; i++)
	{
		Pont >> yPont[i];
		Pont >> xPont[i];
	}

	// your code
	// Definér antal datapunkter
	const int N = 40; // Antal datapunkter
	const int M = 3;  // Antal parametre (a0, a1, a2)

	// Opret designmatrix A og vektor b
	MatDoub A(N, M); // N x M matrix
	VecDoub b(N);	 // N-dimensionel vektor

	// Fyld matricen A og vektoren b
	for (int i = 0; i < N; i++)
	{
		A[i][0] = 1.0;				   // 1. kolonne: konstant led
		A[i][1] = xPont[i];			   // 2. kolonne: x_i
		A[i][2] = xPont[i] * xPont[i]; // 3. kolonne: x_i^2
		b[i] = yPont[i];			   // Højreside-vektor
	}

	// util::print(A);

	// Opret normal matrix A^T * A og A^T * b
	MatDoub ATA(M, M, 0.0);
	VecDoub ATb(M, 0.0);

	// Beregn A^T * A
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < M; j++)
		{
			for (int k = 0; k < N; k++)
			{
				ATA[i][j] += A[k][i] * A[k][j];
			}
		}
	}

	// Beregn A^T * b
	for (int i = 0; i < M; i++)
	{
		for (int k = 0; k < N; k++)
		{
			ATb[i] += A[k][i] * b[k];
		}
	}

	util::print(ATA);

	Cholesky chol(ATA); // Kræver at ATA er positiv definit
	VecDoub a_chol(M);
	chol.solve(ATb, a_chol);

	cout << "Løsning med Cholesky:" << endl;
	for (int i = 0; i < M; i++)
	{
		cout << "a_" << i << " = " << a_chol[i] << endl;
	}



	// FILIP
	cout << "//////////////////////////////////////////" << endl;
	// Definér antal datapunkter
	const int Nf = 82; // Antal datapunkter
	const int Mf = 11;  // Antal parametre (a0, a1, a2)

	// Opret designmatrix A og vektor b
	MatDoub Af(Nf, Mf, 1.0); // N x M matrix
	VecDoub bf(Nf);	 // N-dimensionel vektor

	// Fyld matricen A og vektoren b
	for (int i = 0; i < Nf; i++) {
        double x_power = xFilip[i];  // Starter med x^1
        for (int j = 1; j < Mf; j++) {
            Af[i][j] = x_power;  // x^j
            x_power *= xFilip[i];  // Opdater til x^(j+1)
        }
        bf[i] = yFilip[i];  // Sætter højreside-vektoren
    }

	// util::print(A);

	// Opret normal matrix A^T * A og A^T * b
	MatDoub ATAf(Mf, Mf, 0.0);
	VecDoub ATbf(Mf, 0.0);

	// Beregn A^T * A
	for (int i = 0; i < Mf; i++)
	{
		for (int j = 0; j < Mf; j++)
		{
			for (int k = 0; k < Nf; k++)
			{
				ATAf[i][j] += Af[k][i] * Af[k][j];
			}
		}
	}

	// Beregn A^T * b
	for (int i = 0; i < Mf; i++)
	{
		for (int k = 0; k < Nf; k++)
		{
			ATbf[i] += Af[k][i] * bf[k];
		}
	}

	util::print(ATAf);

	// Cholesky cholf(ATAf); Cholesky fejler fordi den ikke er positiv definit
	LUdcmp lu(ATAf); // Fordi den ikke er positiv definit
	VecDoub a_lu(Mf);
	lu.solve(ATbf, a_lu);

	cout << "Løsning med LU:" << endl;
	for (int i = 0; i < Mf; i++)
	{
		cout << "a_" << i << " = " << a_lu[i] << endl;
	}


	return 0;
}
