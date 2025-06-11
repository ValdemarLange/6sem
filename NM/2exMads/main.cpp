#include "cholesky.h"
#include "ludcmp.h"
#include "nr3.h"
#include "utilities.h"
#include <cmath>
#include <fstream>
#include <iostream>
// your extra include headers

using namespace std;

int main() {
  VecDoub xFilip(82);
  VecDoub yFilip(82);
  ifstream Filip("FilipData.dat");
  for (int i = 0; i < 82; i++) {
    Filip >> yFilip[i];
    Filip >> xFilip[i];
  }

  VecDoub xPont(40);
  VecDoub yPont(40);
  ifstream Pont("PontiusData.dat");
  for (int i = 0; i < 40; i++) {
    Pont >> yPont[i];
    Pont >> xPont[i];
  }

  MatDoub A(40, 3);

  for (int i = 0; i < 40; i++) {
    for (int j = 0; j < 3; j++) {
      A[i][j] = std::pow(xPont[i], j);
    }
  }

  VecDoub b = yPont;

  MatDoub C = util::Transpose(A) * A;
  VecDoub c = util::Transpose(A) * b;

  VecDoub x(3);
  LUdcmp lu(C);
  lu.solve(c, x);
  std::cout << "LU pont" << std::endl;
  util::print(x);

  Cholesky chol(C);
  VecDoub chol_ans(3);
  chol.solve(c, chol_ans);
  std::cout << "Cholesky pont" << std::endl;
  util::print(chol_ans);

  // util::print(chol.el);

  std::cout << "\n";
  std::cout << "Fillip: " << std::endl;

  MatDoub A_fil(82, 11);

  for (int i = 0; i < 82; i++) {
    for (int j = 0; j < 11; j++) {
      A_fil[i][j] = std::pow(xFilip[i], j);
    }
  }

  VecDoub b_fil = yFilip;

  MatDoub C_fil = util::Transpose(A_fil) * A_fil;
  VecDoub c_fil = util::Transpose(A_fil) * b_fil;

  VecDoub x_fil(11);
  LUdcmp lu_fil(C_fil);
  lu_fil.solve(c_fil, x_fil);
  util::print(x_fil);

  Cholesky chol_fil(C_fil);
  VecDoub chol_fil_ans(11);
  chol.solve(c_fil, chol_fil_ans);
  util::print(chol_fil_ans);

  return 0;
}
