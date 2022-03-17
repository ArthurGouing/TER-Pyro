#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <vector>

using namespace std;


int main()
{
  double t=3.;
  std::vector<double> v,rho,x,T;
  double dx=0.0001;
  int n=100;
  double C,rhov,rhop,A,Ta;
  rhov=1500.;
  A=1000.;
  Ta=6000.;
  rhop=1000.; //////t=0; rho=rhop CI
  C=(rhov*A)/(rhov-rhop);
  T.resize(n);
  x.resize(n);
  rho.resize(n);
  for(int i=0;i<n;i++)
  {
    x[i]=i*dx;
    T[i]=200.+80000.*t*(0.01-x[i]);
    rho[i]=rhop*exp(-C*exp(-Ta/T[i])*t)+rhop;
  }

ofstream mon_flux; // Contruit un objet "ofstream"
string name("rho.txt"); // Le nom de mon fichier
mon_flux.open(name,std::ios_base::app | std::ios_base::out); // Ouvre un fichier appelé name_file
if(mon_flux) // Vérifie que le fichier est bien ouvert
{

  for (int i = 0; i <n; i++)
  {
    mon_flux << x[i] << " " << rho[i] << endl;
  }

mon_flux.close();
}
}
