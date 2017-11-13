#include <iostream>
#include <armadillo>
#include <random>
#include <iomanip>

using namespace std;
using namespace arma;
random_device rd;
mt19937 randomEngine(rd());
uniform_real_distribution<double> uniformDist(0.0,1.0);

ofstream outFile;

int n = 2; //size of lattice (nxn)
double J = 1;
double k = 1;
double B,T;


double randomn(){ //returnerer random tall mellom 0 og 1
    return uniformDist(randomEngine);
}

mat lattice(){
    mat lat(n,n,fill::ones);
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if (randomn()<=0.5)
                lat(i,j)*=-1;
        }
    }
    return lat;
}

double energy(mat spin){
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


int main(){

    int M = 1e4; //montecarlosteps
    T=1; //temperatur
    mat spin = lattice();
    double deltaEnergy;
    double Etot  = 0;
    double Etot2 = 0;
    double Magn = 0;
    double Magn2 = 0;
    double count;
    B=1/(T*k);
    count = 0;

 //   string outFileName = "energy_V_" + to_string(n) + "_" + to_string(M) + "_rand_" + to_string(T) + ".txt";
 //   outFile.open(outFileName);


////    string outFileName = "acc_" + to_string(n) + "_" + to_string(M) + "_" + to_string(T) + ".txt";
////    outFile.open(outFileName);

  ///  string outFileName = "accc_" + to_string(n) + "_" + to_string(M) + ".txt";
  ///  outFile.open(outFileName);

        for (int i=0;i<M;i++){
            int x, y;
            pickXY(x, y);
            double Eold = energy(spin);
            spin(x,y) *= -1;
            double Enew = energy(spin);

            deltaEnergy = Enew - Eold;

            if (rejORacc(deltaEnergy) == false) {
                spin(x,y) *= -1;
            }
            else{
                count = count+1;
            }

   ////         outFile << count << endl;
            double E = energy(spin);
            Etot  += E;
            Etot2 += E*E;
            Magn += abs(accu(spin));
            Magn2 += accu(spin)*accu(spin);

        outFile << E << endl;
        }
  //      cout << T << endl;
   ///     outFile << count << "," << T << endl;

 //   outFile.close();

    cout << setprecision(16) << "Energy: " <<Etot/M << endl;
    cout << setprecision(16) << "Magn: " << Magn/M << endl;
    cout << setprecision(16) << "Specific heat: "<< Etot2/(M*n*n) - (Etot/(M*n*n))*(Etot/(M*n*n)) << endl;
    cout << setprecision(16) << "Suceptebility: " << Magn2/(M*n*n) - (Magn/(M*n*n))*(Magn/(M*n*n)) << endl;
   // cout <<"# of changes: " <<count << endl;

    return 0;
}
