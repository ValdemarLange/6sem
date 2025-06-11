#include "nr3.h"
#include "ludcmp.h"
#include "qrdcmp.h"
#include "roots_multidim.h"
#include "utilities.h"
#include <iostream>
#include <math.h>
#include <vector>


VecDoub operator+(VecDoub_I x, VecDoub_I y) {
  VecDoub res(x.size());
  for (int i = 0; i < x.size(); i++) {
    res[i] = x[i] + y[i];
  }

  return res;
}

VecDoub operator-(VecDoub_I x, VecDoub_I y) {
  VecDoub res(x.size());
  for (int i = 0; i < x.size(); i++) {
    res[i] = x[i] - y[i];
  }
  return res;
}

VecDoub operator*(double h, const VecDoub &vec) {
  VecDoub res(vec.size());
  for (int i = 0; i < vec.size(); i++) {
    res[i] = vec[i] * h;
  }

  return res;
}

void printCenteredTitle(std::string title, int width) {
  char fillChar = '-';
  std::cout << std::string(width, fillChar) << std::endl;
  std::cout << std::setw((width + title.length()) / 2) << title << std::endl;
  std::cout << std::string(width, fillChar) << std::endl;
}

void printResults(std::string methodName, std::string x, std::vector<Doub> u,
                  std::vector<Doub> v) {
  std::string title = methodName + " results";
  printCenteredTitle(title, 40);
  std::string title1 = methodName + " results for u(x) with x=" + x;
  std::cout << title1 << std::endl;
  std::cout << std::setw(10) << "n" << std::setw(15) << "A(hi)" << std::setw(15)
            << "A(hi-1)-A(hi)" << std::endl;
  int N = 5;
  for (int i = 0; i < u.size(); i++) {
    if (i == 0) {
      std::cout << std::setw(10) << N << std::setw(15) << u[i] << std::setw(15)
                << "*" << std::endl;
    } else {
      std::cout << std::setw(10) << N << std::setw(15) << u[i] << std::setw(15)
                << u[i - 1] - u[i] << std::endl;
    }
    N *= 2;
  }

  std::cout << "\n";

  std::string title2 = methodName + "results for v(x) with x=" + x;
  std::cout << title2 << std::endl;
  std::cout << std::setw(10) << "n" << std::setw(15) << "A(hi)" << std::setw(15)
            << "A(hi-1)-A(hi)" << std::endl;
  N = 5;
  for (int i = 0; i < v.size(); i++) {
    if (i == 0) {
      std::cout << std::setw(10) << N << std::setw(15) << v[i] << std::setw(15)
                << "*" << std::endl;
    } else {
      std::cout << std::setw(10) << N << std::setw(15) << v[i] << std::setw(15)
                << v[i - 1] - v[i] << std::endl;
    }
    N *= 2;
  }

  std::cout << "\n";
}

VecDoub vecfunc(VecDoub y, double t) {
  VecDoub f(2);

  // y[0] = u
  // y[1] = v
  f[0] = y[0] * y[1];
  f[1] = -(y[0] * y[0]);

  return f;
}

void eulerODE(VecDoub_IO &x, VecDoub (*vecfunc)(VecDoub, double), VecDoub &y0,
              double a, double b, double n) {
  VecDoub ycurr(2);
  ycurr[0] = y0[0];
  ycurr[1] = y0[1];

  VecDoub ynext(2);

  double h = (b - a) / n;

  for (double x = a; x < b; x += h) {
    ynext = ycurr + h * vecfunc(ycurr, x);
    ycurr = ynext;
  }

  x = ycurr;
}

void midpointODE(VecDoub_IO &x, VecDoub (*vecfunc)(VecDoub, double),
                 VecDoub &y0, double a, double b, double n) {

  VecDoub ycurr(y0.size());
  for (int i = 0; i < y0.size(); i++) {
    ycurr[i] = y0[i];
  }

  VecDoub ynext(y0.size());

  double h = (b - a) / n;

  for (double i = a; i < b; i += h) {
    VecDoub k1 = h * vecfunc(ycurr, i);
    VecDoub k2 = h * vecfunc(ycurr + 0.5 * k1, i);
    ynext = ycurr + k2;
    ycurr = ynext;
  }

  x = ycurr;
}

// mintolf er sat til 1e-7 for at kunen køre i newt
void trapezoidODE(VecDoub_IO &x, VecDoub (*vecfunc)(VecDoub, double),
                  VecDoub &y0, double a, double b, double n) {

  VecDoub ycurr(y0.size());
  for (int i = 0; i < y0.size(); i++) {
    ycurr[i] = y0[i];
  }

  VecDoub ynext(y0.size());

  double h = (b - a) / n;

  for (double t = a; t < b; t += h) {
    VecDoub ynext_guess = ycurr + h * vecfunc(ycurr, t);
    Bool check = false;

    auto phi_func = [=](VecDoub ynext_guess) {
      return ynext_guess - ycurr -
             h / 2.0 * (vecfunc(ycurr, t) + vecfunc(ynext_guess, t + h));
    };
    newt(ynext_guess, check, phi_func);
    ycurr = ynext_guess;
  }

  x = ycurr;
}

void print(VecDoub_I &vec, const string &text = "") {
  cout << text << "  Vector " << vec.size() << "D:" << endl;
  for (int m = 0; m < vec.size(); m++) {
      cout << setw(15) << vec[m];
  }
  cout << endl;
}

void fourthOrderRungeKutta(VecDoub_IO &x, VecDoub (*vecfunc)(VecDoub, double),
                           VecDoub &y0, double a, double b, double n) {
  VecDoub ycurr(y0.size());
  for (int i = 0; i < y0.size(); i++) {
    ycurr[i] = y0[i];
  }

  VecDoub ynext(y0.size());
  double h = (b - a) / n;

  for (double t = a; t < b; t += h) {
    VecDoub k1 = h * vecfunc(ycurr, t);
    VecDoub k2 = h * vecfunc(ycurr + 0.5 * k1, t);
    VecDoub k3 = h * vecfunc(ycurr + 0.5 * k2, t);
    VecDoub k4 = h * vecfunc(ycurr + k3, t);
    ynext = ycurr + 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    ycurr = ynext;
  }
  print(ynext, "ynext");
  x = ycurr;
}

void printEuler(VecDoub &y0) {
  VecDoub_IO x(2);

  std::vector<Doub> u;
  std::vector<Doub> v;

  int N = 5;
  for (int i = 0; i < 10; i++) {
    eulerODE(x, vecfunc, y0, 0, 10, N);
    u.push_back(x[0]);
    v.push_back(x[1]);
    N *= 2;
  }

  printResults("Euler", "10", u, v);
}

void printMid(VecDoub &y0) {
  VecDoub_IO x(2);

  std::vector<Doub> u;
  std::vector<Doub> v;

  int N = 5;
  for (int i = 0; i < 10; i++) {
    midpointODE(x, vecfunc, y0, 0, 10, N);
    u.push_back(x[0]);
    v.push_back(x[1]);
    N *= 2;
  }

  printResults("Midpoint", "10", u, v);
}

void printTrapezoidal(VecDoub &y0) {
  VecDoub_IO x(2);
  std::vector<Doub> u, v;
  int N = 5;
  for (int i = 0; i < 10; i++) {
    trapezoidODE(x, vecfunc, y0, 0, 10, N);
    u.push_back(x[0]);
    v.push_back(x[1]);
    N *= 2;
  }
  printResults("Trapezoidal", "10", u, v);
}

void printFourthRunge(VecDoub &y0) {
  VecDoub_IO x(2);
  std::vector<Doub> u, v;
  int N = 5;
  for (int i = 0; i < 10; i++) {
    fourthOrderRungeKutta(x, vecfunc, y0, 0, 10, N);
    u.push_back(x[0]);
    v.push_back(x[1]);
    N *= 2;
  }
  printResults("Fourth order runge kutta", "10", u, v);
}

int main() {
  VecDoub y0(2);
  y0[0] = 1;
  y0[1] = 1;

  // printEuler(y0);
  // printMid(y0);
  // printTrapezoidal(y0);
  printFourthRunge(y0);

  return 1;
}
