#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"


using namespace Eigen;
using namespace std;

int main(int argc, char** argv) // ./laplacian dataSmallCase.toml -> argc=2 et argv[0]="./laplacian" argv[1]="dataSmallCase.toml"
{
  if (argc < 2)
  {
  cout << "Please, enter the name of your data file." << endl; abort();
  }
  // ----------------------- Fichier de données --------------------------------
  DataFile* data_file = new DataFile(argv[1]);

  // ---------------------------- Résolution  ----------------------------------
  double nb_iterations=data_file->Get__nb_iterations();
  // Pointeur contenant toutes les fonctions utiles
  Function* function = new Function(data_file);
  // Pointeur vers la classe MeshAdapt (adaptation de maillage) //!!!!!
  cout << "constructeur"<<endl;
  Mesh_Adapt* mesh_adapt = new Mesh_Adapt(data_file); //!!!!
  // Pointeur vers la classe FiniteVolume (discrétisation en espace)
  FiniteVolume* fin_vol = new FiniteVolume(function, data_file, mesh_adapt);
  // Pointeur vers la classe TimeScheme (discrétisation en temps=)
  TimeScheme* time_scheme = NULL;
  time_scheme = new ImplicitEulerScheme(data_file, fin_vol);

  Eigen::VectorXd sol;
  double tn, dt=data_file->Get_dt();
  string scenario=data_file->Get_scenario();

  cout << "-------------------------------------------------" << endl;
  cout << "Search u such that : " << endl;
  cout << "dt u + div (v u) - div(mu grad(u)) = f" << endl;
  cout << "-------------------------------------------------" << endl;

  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();

  string n_file0 ="results/results0.dat", n_file1 ="results/results1.dat", n_file2 ="results/results2.dat", n_file3 ="results/results3.dat", n_file4 ="results/results4.dat";
  ofstream solu0, solu1, solu2, solu3, solu4;
  solu0.open(n_file0, ios::out);
  solu1.open(n_file1, ios::out);
  solu2.open(n_file2, ios::out);
  solu3.open(n_file3, ios::out);
  solu4.open(n_file4, ios::out);
  int Nx=data_file->Get_Nx(), Ny=data_file->Get_Ny();
  double dy = data_file->Get_dy();
  int Ny0(Nx*(Ny-int(ceil(0./dy)))-1), Ny1(Nx*(Ny-int(ceil(0.001/dy)))-1), Ny2(Nx*(Ny-int(ceil(0.002/dy)))-1), Ny3(Nx*(Ny-int(ceil(0.003/dy)))-1), Ny4(Nx*(Ny-int(ceil(0.004/dy)))-1);
  cout << Ny0 << " " << Ny1 << " " << Ny2 << " " << Ny3 << " " << Ny4 << endl;

  cout << "Save initial condition " << endl;
  time_scheme->SaveSol(time_scheme->GetSolution(),"ImpliciteScheme", 0); //peut être que ca va déconner avec le fait de pas utiliser une varibale
  cout << "Time Loop" << endl;
  sol=time_scheme->GetSolution();
  solu0 << 0. << " " << sol(Ny0) << endl;
  solu1 << 0. << " " << sol(Ny1) << endl;
  solu2 << 0. << " " << sol(Ny2) << endl;
  solu3 << 0. << " " << sol(Ny3) << endl;
  solu4 << 0. << " " << sol(Ny4) << endl;


  for (int n = 1; n <= nb_iterations; n++) // Boucle en temps
  {
    cout << "Iteration " << n << endl;
    time_scheme->Advance();
    //time_scheme->SaveSol(time_scheme->GetSolution(),"ImpliciteScheme", n);
    sol=time_scheme->GetSolution();

    //Savoir comment faire ????
    if (scenario=="adapt")
    {
    mesh_adapt->Update(sol);
    }

    tn=n*dt;
    solu0 << tn << " " << sol(Ny0) << endl;
    solu1 << tn << " " << sol(Ny1) << endl;
    solu2 << tn << " " << sol(Ny2) << endl;
    solu3 << tn << " " << sol(Ny3) << endl;
    solu4 << tn << " " << sol(Ny4) << endl;
  }
  time_scheme->SaveSol(time_scheme->GetSolution(),"ImpliciteScheme", nb_iterations);
  // ---------------------------------------------------------------------------

  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::seconds>(finish-start).count();
  // Affichage du résultat
  cout << "Cela a pris "<< t << " seconds" << endl;

  delete time_scheme;
  delete fin_vol;
  delete data_file;
  delete function;
  delete mesh_adapt;
  return 0;
}
