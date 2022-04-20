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
  Mesh_Adapt* mesh_adapt = new Mesh_Adapt(data_file); //!!!!
  // Pointeur vers la classe FiniteVolume (discrétisation en espace)
  FiniteVolume* fin_vol = new FiniteVolume(function, data_file, mesh_adapt);
  // Pointeur vers la classe TimeScheme (discrétisation en temps=)
  TimeScheme* time_scheme = NULL;
  time_scheme = new ImplicitEulerScheme(data_file, fin_vol);

  Solution Solution_tn(data_file);
  Eigen::VectorXd temp, rho;
  double tn=data_file->Get_t0(), dt=data_file->Get_dt();
  string maillage=data_file->Get_scenario(), scenario=data_file->Get_results();

  cout << "                                                  " << ::endl;
  cout << "--------------------------------------------------" << endl;
  cout << "------------ Search T and p such that : ----------" << endl;
  cout << " dt(p*Cp*(T-T0)) = div(lambda*grad(T)) - Lm*dt(p) " << endl;
  cout << "    dt(p) = -pv*A*((p-pp)/(pv-pp))*exp(-Ta/T)     " << endl;
  cout << "------------ Case : "<< scenario << " " << maillage << " -------------" << endl;
  cout << "--------------------------------------------------" << endl;


  // Démarrage du chrono
  auto start = chrono::high_resolution_clock::now();


  //Ouverture des fichiers solutions à plusieurs distances de la paroie
  vector<ofstream> templist;
  vector<ofstream> rholist;
  vector<ofstream> tempylist;
  templist.resize(5);
  rholist.resize(5);
  tempylist.resize(11);
  for (int i = 0; i <= 4 ; i++)
  {
    string n_filetemp = data_file->Get_results() + "Temperatures_" + to_string(i)+ "mm.dat";
    templist[i].open(n_filetemp, ios::out);
    if (data_file->Get_Aref()!=0)
    {
      string n_filerho  = data_file->Get_results() + "Rho_" + to_string(i)+ "mm.dat";
      rholist[i].open(n_filerho, ios::out);
    }
  }


  //Température suivant y à t0
  string n_filetemp = data_file->Get_results() + "Temperatures_suivant_y_a_" + to_string(data_file->Get_t0())+ "_secondes.dat";
  tempylist[0].open(n_filetemp, ios::out);


  //Obtention du numéro de la ligne corespondant à la distance
  vector<int> Nylist(5);
  Nylist={mesh_adapt->cellule(0.),mesh_adapt->cellule(0.001),mesh_adapt->cellule(0.002),mesh_adapt->cellule(0.003),mesh_adapt->cellule(0.004)};
  cout << "                                                  " << ::endl;
  cout << "--------------------------------------------------" << endl;
  cout << "   Numéro de la ligne correpondant à la distance  " << endl;
  cout << "A 0 mm : Ny0 = " << Nylist[0] << "        " << "A 1 mm : Ny1 = " << Nylist[1] << endl;
  cout << "A 2 mm : Ny2 = " << Nylist[2] << "        " << "A 3 mm : Ny3 = " << Nylist[3] << endl;
  cout << "A 4 mm : Ny4 = " << Nylist[4] << endl;
  cout << "!!! Si vous obtenez des nombres négatifs : ymax est trop petit par rapport à la distance à la paroi où on cherche la solution !!!"<<endl;
  cout << "--------------------------------------------------" << endl;

  //Sauvegarde de la solution initiale
  cout << "                                                  " << ::endl;
  cout << "--------------------------------------------------" << endl;
  cout << "------------ Save initial condition --------------" << endl;
  time_scheme->SaveSol(time_scheme->Get_Solution(),"ImpliciteScheme", 0); //peut être que ca va déconner avec le fait de pas utiliser une varibale
  temp=time_scheme->Get_Solution().Get_T();
  rho=time_scheme->Get_Solution().Get_rho();
  Solution_tn=time_scheme->Get_Solution();
  for (int i = 0; i <= 4 ; i++)
  {
    templist[i] << 0. << " " << temp(Nylist[i]) << endl;
    if (data_file->Get_Aref()!=0.)
    {
      rholist[i]  << 0. << " " << rho(Nylist[i])  << endl;
    }
  }
  cout << "-------------------------------------------------" << endl;


  //Température suivant y à t0
  double dist = (mesh_adapt->Get_Dy())(0)/2;
  tempylist[0] << dist << " " << temp(0) << endl;
  int Ny = data_file->Get_Ny();
  int Nx = data_file->Get_Nx();
  for (int i=1; i<Ny; i++)
  {
    dist += ((mesh_adapt->Get_Dy())(i-1)+(mesh_adapt->Get_Dy())(i))/2;
    tempylist[0] << dist << " " << temp(i*Nx) << endl;
  }
  cout << "-------------------------------------------------" << endl;


  //Itérations
  cout << "                                                  " << ::endl;
  cout << "--------------------------------------------------" << endl;
  cout << "------------      Time Loop         --------------" << endl;
  double dtgraph = (data_file->Get_tfinal()/10); //10 temps differents pour la température suivant y ??
  double cond(dtgraph);
  int compteur(0);
  for (int n = 1; n <= nb_iterations; n++) // Boucle en temps
  {
    if (maillage == "noadapt")
    {
      tn=n*dt;
      cout << "Iteration : " << n << " Temps : " << tn << " s" << endl;
      time_scheme->Advance();
    }
    else if (maillage == "adapt")
    {
      tn=n*dt;
      cout << "Iteration : " << n << " Temps : " << tn << " s" << endl;
      mesh_adapt->Update_Dyold();
      time_scheme->Update_Told_rhoold();                                        // Il faut faire une sauvegarde de Told=Tn et rhoold=rhon (CI)
      time_scheme->Advance();                                                   // On effectue une premiere initialisation de Tn+1 et rhon+1

      //Adaptation de maillage
      mesh_adapt->Update_Dyprevious();                                          // On sauvegarde Dyk(n) pour le calcul de la norme : |Dyk(n)-Dyk+1(n)|
      cout << "1" << endl;
      mesh_adapt->Update(time_scheme->Get_Solution());                          // Calcul du premier maillage adapté
      cout << "2" << endl;
      mesh_adapt->Update_Dystar_vitesse();                                      // Calcul de Dystar et de la vitesse d'advection du maillage
      cout << "3" << endl;
      time_scheme->Advance_ALE();                                               // Calcul de Tn+1 sur nouveau maillage
      mesh_adapt->Update_Dyold();                                             ////Peut etre a retiré
    //   //Pour éviter absolument de faire une initialisation pour intialiser Normlinf en dehors de la boucle on peut faire poser une variable norme et la mettre à eps+1 au début
    //   cout << mesh_adapt->NormLinf()<< " " <<data_file->Get_epsilon_adapt() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl;
    //   while (mesh_adapt->NormLinf() > data_file->Get_epsilon_adapt())           // Boucle pour rafinner le maillage
    //   {
    //     mesh_adapt->Update_Dyprevious();                                        // On sauvegarde Dy(n) pour le calcul de la norme : |Dyk(n)-Dyk+1(n)|
    //     mesh_adapt->Update(time_scheme->Get_Solution());                        // Calcul de Dyk+1(n)
    //     mesh_adapt->Update_Dystar_vitesse();                                    // Calcul de Dystar et de la vitesse d'advection du maillage
    //     time_scheme->Advance_ALE();                                             // Calcul de Tn+1 sur nouveau maillage
    //     cout << mesh_adapt->NormLinf()<< " " << data_file->Get_epsilon_adapt() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<endl;
    //   }
    //   cout << "après boucle adapt maillage" << endl;
    //   mesh_adapt->Update_Dyold();                                               // Le dernier maillage de la boucle devient le Dyold pour la CI pour la nouvelle itération (Tn+2)
    }


    // Ecriture
    time_scheme->SaveSol(time_scheme->Get_Solution(),"ImpliciteScheme", n);
    temp = time_scheme->Get_Solution().Get_T();                                 // Sauvegarde des vecteurs pour écriture
    rho  = time_scheme->Get_Solution().Get_rho();

    for (int i = 0; i <= 4 ; i++)
    {
      templist[i] << tn << " " << temp(Nylist[i]) << endl;
      if (data_file->Get_Aref()!=0.)
      {
        rholist[i]  << tn << " " << rho(Nylist[i])  << endl;
      }
    }

    //Temperature suivant y
    if (abs(tn-cond) <= dt/2)
    {
      cond+=dtgraph;
      compteur+=1;
      dist = (mesh_adapt->Get_Dy())(0)/2;
      n_filetemp = data_file->Get_results() + "Temperatures_suivant_y_a_" + to_string(tn)+ "_secondes.dat";
      tempylist[compteur].open(n_filetemp, ios::out);
      tempylist[compteur] << dist << " " << temp(0) << endl;
      for (int j=1; j<Ny; j++)
      {
        dist += ((mesh_adapt->Get_Dy())(j-1)+(mesh_adapt->Get_Dy())(j))/2;
        tempylist[compteur] << dist << " " << temp(j*Nx) << endl;
      }
    }

  }
  cout << "--------------------------------------------------" << endl;

  // Fin du chrono
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::seconds>(finish-start).count();
  // Affichage du résultat
  cout << "Cela a pris "<< t << " secondes" << endl;

  delete time_scheme;
  delete fin_vol;
  delete data_file;
  delete function;
  delete mesh_adapt;
  return 0;
}
