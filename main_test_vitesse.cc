#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <vector>

using namespace std;

int main()
{
  std::vector<double> v, Y,Y_1;
  int n=10;
  Y.resize(n);
  Y_1.resize(n);
  for (int i=0;i<n;i++)
  {
    Y[i]=1+i;
  }
  Y_1=Y;

  for (int i=0;i<n;i++)
  {
    Y[i]=10;
  }
  v.resize(Y_1.size());
  for(int i=0; i<v.size();i++)
  {
    v[i]=(Y[i+1]-Y[i])-(Y_1[i+1]-Y_1[i]);
  }


  ofstream mon_flux; // Contruit un objet "ofstream"
  string name("dina.txt"); // Le nom de mon fichier
  mon_flux.open(name,std::ios_base::app | std::ios_base::out); // Ouvre un fichier appelé name_file
  if(mon_flux) // Vérifie que le fichier est bien ouvert
  {

    for (int i = 0; i <n; i++)
    {
      mon_flux << Y[i] << " " << Y_1[i] << endl;
    }

  mon_flux.close();
  }





}
