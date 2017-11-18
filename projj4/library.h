#ifndef LIBRARY_H
#define LIBRARY_H
#include <armadillo>

/*
void pickXY(int& x, int& y, int n);
bool rejORacc(double deltaE, int B);
void lattice(arma::mat& spin, int n);
double energy(arma::mat& spin, int n, int J,double& Magnum);
*/


void pickXY(int& x, int& y);
bool rejORacc(double deltaE);
arma::mat lattice();
double energy(arma::mat& spin);
void lattice(arma::mat& spin, int n);
#endif // LIBRARY_H
