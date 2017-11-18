#include <iostream>
#include <armadillo>
#include <random>
#include <iomanip>
#include "library.h"
using namespace std;
using namespace arma;
ofstream outFile;

int n = 20; //size of lattice (nxn)
double J = 1;
double k = 1;
double B,T;

int main(){

    int M = 1e5*n*n; //montecarlosteps
    T=2.4; //temperatur
    mat spin = ones(n,n);
    lattice(spin,n);
    double deltaEnergy;
    double Etot  = 0;
    double Etot2 = 0;
    double Magn = 0;
    double Magn2 = 0;
    double count;
    B=1/(T*k);
    count = 0;

    string outFileName = "cummuEran_" + to_string(n) + "_" + to_string(T) + ".txt";
    outFile.open(outFileName);

    double E = energy(spin);
        for (int i=0;i<M;i++){
            int x, y;
            pickXY(x, y);

            int xv,xh,yn,yu;
            xh=(x+1)%n; //høyre
            yn=(y+1)%n; //ned
            xv=(x-1+n)%n; //venstre
            yu=(y-1+n)%n; //opp


            deltaEnergy = 2*J*(spin(x,y)*(spin(xh,y)+spin(x,yn)+spin(xv,y)+spin(x,yu)));

            if (rejORacc(deltaEnergy) == true) {
                spin(x,y) *= -1;
                E+=deltaEnergy;
                count = count+1;
                Magn +=2*spin(x,y);

            }
            Etot  += E;
            Etot2 += E*E;
            Magn2 += Magn*Magn;

            if (i%(n*n)==0){
                outFile << Etot << ", " << i << endl;
            }

        }

    outFile.close();

    cout << setprecision(16) << "Energy: " << Etot/M/n/n << endl;
    cout << setprecision(16) << "Magn: " << Magn/M << endl;
    cout << setprecision(16) << "Specific heat: "<< Etot2/(M*n*n) - (Etot/(M*n*n))*(Etot/(M*n*n)) << endl;
    cout << setprecision(16) << "Suceptebility: " << Magn2/(M*n*n) - (Magn/(M*n*n))*(Magn/(M*n*n)) << endl;
    cout <<"# of changes: " <<count << endl;

    return 0;
}
