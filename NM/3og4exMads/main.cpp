#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include "nr3.h"
#include "ludcmp.h"
#include "utilities.h"
#include "svd.h"
// your extra include headers

using namespace std;

VecDoub print_error(SVD SVD_A, std::string name) {
  MatDoub V = SVD_A.v;
  VecDoub W = SVD_A.w;

  VecDoub sigma(W.size());
  double sum;
  for (int j = 0; j < W.size(); j++) {
    sum = 0;
    for (int i = 0; i < W.size(); i++) {
      if (W[j] > pow(10, -15)) {
        sum += pow(V[j][i] / W[i], 2);
      }
    }
    sigma[j] = sqrt(sum);
  }
  std::string print_name =
      name + " standard deviation of estimated model parameters a:";
  util::print(sigma, print_name);
  return sigma;
}

double vectorLength(VecDoub v) {
  double sum = 0;
  for (int i = 0; i < v.size(); i++) {
    sum += pow(v[i], 2);
  }
  return sqrt(sum);
}

double print_residual_error(MatDoub A, VecDoub x, VecDoub b, std::string name) {
  VecDoub dif = A * x - b;
  Doub n = A.ncols();
  Doub m = A.nrows();
  Doub randomFit = sqrt((m - n) / m);
  double residual_error = vectorLength(dif) / vectorLength(b);
  std::cout << "Standard residual error of " << name << ": " << residual_error
            << " - Compared to random fit error: " << randomFit << std::endl;
  return residual_error;
}

void vecAbs(VecDoub &v) {
  for (int i = 0; i < v.size(); i++) {
    v[i] = abs(v[i]);
  }
}

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

  std::cout << "Pont: " << std::endl;

  VecDoub x(3);
  SVD sv(A);
  sv.solve(b, x);
  util::print(x);

  VecDoub sigma_pont = print_error(sv, "Pont");
  double residual_pont = print_residual_error(A, x, b, "Pont");

  std::cout << "RANK: " << sv.rank() << std::endl;

  std::cout << "\n";
  std::cout << "Pontius with stated std accounted for" << std::endl;
  MatDoub A_pont_std(40, 3);
  VecDoub b_pont_std(40);

  VecDoub std_pont(40);
  std_pont = A * x - b;

  double curr_sigma = residual_pont;

  for (int i = 0; i < 40; i++) {
    for (int j = 0; j < 3; j++) {
      // curr_sigma = max(delta, abs(std_pont[i]));
      A_pont_std[i][j] = A[i][j] / curr_sigma;
    }
    b_pont_std[i] = b[i] / curr_sigma;
  }

  SVD sv_pont_std(A_pont_std);
  VecDoub x_pont_std(3);
  sv_pont_std.solve(b_pont_std, x_pont_std);

  // TODO: check if the correction works, seems to be good :)
  util::print(x_pont_std, "Estimated model parameters with std accounted for");
  residual_pont =
      print_residual_error(A_pont_std, x_pont_std, b_pont_std, "Pontuis");
  print_error(sv_pont_std, "Pontius");

  std::cout << "\n";
  std::cout << "Fillip: " << std::endl;

  MatDoub A_fil(82, 11);

  for (int i = 0; i < 82; i++) {
    for (int j = 0; j < 11; j++) {
      A_fil[i][j] = std::pow(xFilip[i], j);
    }
  }

  VecDoub b_fil = yFilip;
  // Set threshold as low as numerical limit
  Doub threshold = std::numeric_limits<float>::denorm_min();
  VecDoub x_fil(11);
  SVD sFil(A_fil);
  sFil.solve(b_fil, x_fil, threshold);
  util::print(x_fil);

  VecDoub sigma_fil = print_error(sFil, "Filip");
  double residual_fil = print_residual_error(A_fil, x_fil, b_fil, "Pont");

  std::cout << "\n";
  std::cout << "Filip with stated std accounted for" << std::endl;
  MatDoub A_fil_std(82, 11);
  VecDoub b_fil_std(82);

  curr_sigma = residual_fil;

  for (int i = 0; i < 82; i++) {
    for (int j = 0; j < 11; j++) {
      A_fil_std[i][j] = A_fil[i][j] / curr_sigma;
    }
    b_fil_std[i] = b_fil[i] / curr_sigma;
  }

  SVD sv_fil_std(A_fil_std);
  VecDoub x_fil_std(11);
  sv_fil_std.solve(b_fil_std, x_fil_std);
  util::print(x_fil_std, "Estimated model parameters with std accounted for");
  print_residual_error(A_fil_std, x_fil_std, b_fil_std, "Filip");
  print_error(sv_fil_std, "Filip");

  // std::cout << "Fillip range: ";
  // MatDoub range = sFil.range(1e-15);

  return 0;
}
