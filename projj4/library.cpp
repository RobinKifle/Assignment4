#include <iostream>
#include "library.h"
#include <random>
#include <iomanip>

using namespace std;
using namespace arma;
random_device rd;
mt19937 randomEngine(rd());
uniform_real_distribution<double> uniformDist(0.0,1.0);

extern int n;
extern double J;
extern double k;
extern double B,T;


double randomn(){ //returnerer random tall mellom 0 og 1
    return uniformDist(randomEngine);
}

void lattice(mat& spin, int n){ //randomizes matrix
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if (randomn()<=0.5)
                spin(i,j)*=-1;
        }
    }
}

double energy(mat& spin){
    double xh, yn, deltaE, Esum;
    Esum = 0;

    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){

            xh=(i+1+n)%n; //høyre
            yn=(j+1+n)%n; //ned

            deltaE = -J*(spin(i,j)*(spin(xh,j)+spin(i,yn)));

            Esum+=deltaE;
        }
    }
    return Esum;
}

void pickXY(int& x, int& y){ //returnerer random point(x,y)
    x = (int) (randomn()*n);
    y = (int) (randomn()*n);
}



bool rejORacc(double deltaE){
    if (exp(-B*deltaE)>=randomn())
        return true;
    else
        return false;
}
