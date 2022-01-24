#ifndef _MESHADAPT_CPP


#include "MeshAdapt.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen"

using namespace Eigen;
using namespace std;


Mesh_Adapt::Mesh_Adapt(DataFile* data_file) :
_df(data_file)
{
  _Dy.resize(_df->Get_Ny()); //taille de Dy=nombre de cases verticales
  _Y.resize(_df->Get_Ny()+1);
  _Y(0)=0.;
  for (int i=0; i<_Dy.size(); i++)
  {
    _Dy(i) =_df->Get_dy();
    _Y(i+1) = _Y(i)+_df->Get_dy();
  }

  ////////////////////////////Nouveau pour le test
  int i=0;
  ifstream fichier("Eulerexplicite_ci_1_cl_1_L_0._tmax_4.0_imax_1000.dat", ios::in);  //Ouverture d'un fichier en lecture
  if(fichier)
  {
    string ligne; //Une variable pour stocker les lignes lues
    double number1, number2;

    _Y.resize(1001);
    _T.resize(1001);
    _Dy.resize(1000);
    while (fichier >> number1 >> number2)
    {
      _Y(i)=number1;
      _T(i)=number2;
      cout << _Y(i) << " " << _T(i) << endl;
      i=i+1;
    }
    // while(getline(fichier, ligne)) //Tant qu'on n'est pas à la fin, on lit
    // {
    //   double a, b;
    //    cout << ligne << endl;
    //    fichier >> a >> b;
    //    cout << a << b << endl;
    //    // fichier >>_Y(i) >> _T(i);
    //    // cout << _Y(i) << _T(i);
    //    //Et on l'affiche dans la console
    //    //Ou alors on fait quelque chose avec cette ligne
    //    //À vous de voir
    //  }

  }
  else
  {
    cout << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
  }

}

void Mesh_Adapt::Update(VectorXd T)
{
  //calcul u
  double Nx = _df->Get_Nx();
  double Ny = 1000;//= _df->Get_Ny();
  double Ly= Ny*_df->Get_dy();
  VectorXd U2(_Dy.size()+1);
  _T= T;

  U2=Derive_y_2(_T); // dérivée seconde selon y en x = dx, aux noeud du maillaage // Attention aux bord
  for (int i=0 ; i<U2.size(); i++)
  {
    cout << U2(i) << " " <<  pow(40.,2)*exp(40.*_Y(i)) << endl;
  }
  //calcul de k
  VectorXd metric(Ny+1); // pour pas faire 2 fois le meme calcul de metric
  VectorXd K(Ny); //vecteur inconnu x(1 à Ny-1)

  for (int i=0 ; i<metric.size(); i++)  //Calcul de la métrique
  {
    double minimum = 1.;
    metric(i)= sqrt( max(U2(i),0.1));
  }

  save_vector(metric,_Y, "test_Metric");

  for (int i=0 ; i<K.size(); i++)  //Calcul de K
  {
    K(i)=(metric(i)+metric(i+1))/2.;
  }

  save_vector(K,_Dy, "test_K");
  //on construit M et b puis on résout M.x = b
  SparseMatrix<double> M(_Y.size()-2, _Y.size()-2);
  VectorXd b(_Y.size()-2), Yprime;
  vector<Triplet<double>> triplets;
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;

  save_vector_mesh(_Y, "test_y0");
  save_vector(U2,_Y, "test_U2");
  b.setZero();
  b(b.size()-1)=K(K.size()-1)*Ly;

  triplets.push_back({0,0,K(1)+K(0)});
  triplets.push_back({0,1,-K(1)});


  for (int i=1; i<M.rows()-1;i++) //On remplit M
  {
    triplets.push_back({i,i-1,-K(i)});
    triplets.push_back({i,i,K(i)+K(i+1)});
    triplets.push_back({i,i+1,-K(i+1)});
  }

  triplets.push_back({M.rows()-1,M.rows()-2,-K(M.rows()-1)});
  triplets.push_back({M.rows()-1,M.rows()-1,K(M.rows()-1)+K(M.rows())});

  M.setFromTriplets(triplets.begin(), triplets.end());

  solver.analyzePattern(M);
  solver.factorize(M);

  Yprime = solver.solve(b);

  for (int i=1; i<_Y.size()-1;i++) //On remplit M
  {
    _Y(i)=Yprime(i-1);
  }

  // On détermine les dy a partir de des Xi
  //_Dy(0)=_Y(0);
  //_Dy(_Dy.size()-1)=Ly-_Y(_Y.size()-2);//Ly-_Y(_Y.size()-1);
  //for (int i=1;i<_Dy.size()-1;i++)
  for (int i=0;i<_Dy.size();i++)
  {
    _Dy(i)= _Y(i+1)-_Y(i);
  }

  save_vector_mesh(_Y, "test_Y_update");

  cout <<"marche !!!"<<endl;
}


VectorXd Mesh_Adapt:: Derive_y_2(VectorXd T)
{
  double Ny = 1000;//_df->Get_Ny();
  double Nx = _df->Get_Nx();
  double dy= _df->Get_dy();
  VectorXd derive2(Ny+1), derive1(Ny+1);
  VectorXd U(Ny+1);

  //Calcul de U
  U(0)=exp(0);
  for (int i=1;i<=Ny-1; i++)
  {
    U(i)=exp(40.*_Y(i)); //T(i-1)+T(i)
  }
  U(Ny)=exp(40.*_Y(Ny));

    //Calcul de U
    // U(0)=T(0);
    // for (int i=1;i<=Ny-1; i++)
    // {
    //   U(i)=(T((i-1)*Nx)+T(i*Nx))/2.;
    // }
    // U(Ny)=T((Ny-1)*Nx);

  //Calcul de U1
  // derive1(0)=(U(1)-U(0))/(_Y(1)-_Y(0));
  // for (int i=1;i<=Ny-1; i++)
  // {
  //   derive1(i)=(U(i+1)-U(i-1))/(_Y(i+1)-_Y(i-1));
  // }
  // derive1(Ny)=(U(Ny)-U(Ny-1))/(_Y(Ny)-_Y(Ny-1));

  //Calcul de U2
  for (int i=1;i<=Ny-1; i++) //Cas général (1 à Ny-1)
  {
    derive2(i) = 2*(_Dy(i-1)*U(i+1) + _Dy(i)*U(i-1) - (_Dy(i-1)+_Dy(i))*U(i)) / (_Dy(i)*_Dy(i)*_Dy(i-1) + _Dy(i-1)*_Dy(i-1)*_Dy(i));
  }
  derive2(0) = derive2(1);
  derive2(Ny) = derive2(Ny-1);


  // //Calcul de U1
  // derive1(0)=(U(1)-U(0))/(_Y(1)-_Y(0));
  // cout << "dérivée1" << derive1(0) << endl;
  // for (int i=1;i<=Ny-1; i++)
  // {
  //   derive1(i)=(U(i+1)-U(i-1))/(_Y(i+1)-_Y(i-1));
  // }
  // derive1(Ny)=(U(Ny)-U(Ny-1))/(_Y(Ny)-_Y(Ny-1));
  // cout << "dérivée1Ny" << derive1(Ny) << endl;
  //
  // //Calcul de U2
  // cout <<  (2./pow(_Y(1)-_Y(0),2)) << " " << (U(1)-U(0)) << " " <<(_Y(1)-_Y(0))*derive1(0) << endl;
  // derive2(0) = (2./pow(_Y(1)-_Y(0),2))*(U(1)-U(0)-(_Y(1)-_Y(0))*derive1(0));
  // cout << "dérivée2" << " " << derive2(0) << endl;
  // derive2(0) = (derive1(1)-derive1(0))/(_Y(1)-_Y(0));;
  // cout << "dérivée2" << " " << derive2(0) << endl;
  // for (int i=1;i<=Ny-1; i++) //Cas général (1 à Ny-1)
  // {
  //   double U2g = (2./pow(_Y(i-1)-_Y(i),2))*(U(i-1)-U(i)-(_Y(i-1)-_Y(i))*derive1(i));
  //   double U2d = (2./pow(_Y(i+1)-_Y(i),2))*(U(i+1)-U(i)-(_Y(i+1)-_Y(i))*derive1(i));
  //
  //   derive2(i) = (U2d*(_Y(i+1)-_Y(i))-U2g*(_Y(i-1)-_Y(i)))/((_Y(i+1)-_Y(i))-(_Y(i-1)-_Y(i)));
  // }
  // cout << "what" << U(Ny-1) << U(Ny) << endl;
  // derive2(Ny) = (2./pow(_Y(Ny-1)-_Y(Ny),2))*(U(Ny-1)-U(Ny)-(_Y(Ny-1)-_Y(Ny))*derive1(Ny));


  return derive2;
}


// VectorXd Mesh_Adapt:: Derive_y_2(VectorXd T)
// {
//   double Ny = _df->Get_Ny();
//   double Nx = _df->Get_Nx();
//   VectorXd derive2(Ny+1), derive1(Ny+1);
//   VectorXd U(Ny+1);
//
//   //Calcul de U
//   U(0)=T(0);
//   for (int i=1;i<=Ny-1; i++)
//   {
//     U(i)=(T((i-1)*Nx)+T(i*Nx))/2.;
//   }
//   U(Ny)=T((Ny-1)*Nx);
//
//   //Calcul de U1
//   derive1(0)=(U(1)-U(0))/(_Y(1)-_Y(0));
//   for (int i=1;i<=Ny-1; i++)
//   {
//     U(i)=(T(i-1)+T(i))/2.;
//   }
//   derive1(Ny)=(U(Ny)-U(Ny-1))/(_Y(Ny)-_Y(Ny-1));
//
//   //Calcul de U2
//   derive2(0) = (2./pow(_Y(1)-_Y(0),2))*(U(1)-U(0)-(_Y(1)-_Y(0))*derive1(0));
//   for (int i=1;i<=Ny-2; i++) //Cas général (1 à Ny-1)
//   {
//     double U2g = (2./pow(_Y(i-1)-_Y(i),2))*(U(i-1)-U(i)-(_Y(i-1)-_Y(i))*derive1(i));
//     double U2d = (2./pow(_Y(i+1)-_Y(i),2))*(U(i+1)-U(i)-(_Y(i+1)-_Y(i))*derive1(i));
//
//     derive2(i) = (U2d*(_Y(i+1)-_Y(i))-U2g*(_Y(i-1)-_Y(i)))/((_Y(i+1)-_Y(i))-(_Y(i-1)-_Y(i)));
//   }
//   derive2(Ny) = (2./pow(_Y(Ny-1)-_Y(Ny),2))*(U(Ny-1)-U(Ny)-(_Y(Ny-1)-_Y(Ny))*derive1(Ny));
//
//   return derive2;
// }



void Mesh_Adapt::save_vector_mesh(Eigen::VectorXd Y, std::string a) // pour le mesh
{
  ofstream flux;
  flux.open(a);
  for (int i=0; i<Y.size();i++)
  {
    flux << Y(i) << " " << 0 << endl;
  }
  flux.close();
}

void Mesh_Adapt::save_vector(Eigen::VectorXd U, Eigen::VectorXd Y, std::string a) // pour voir U2
{
  ofstream flux;
  flux.open(a);
  for (int i=0; i<Y.size();i++)
  {
    flux << Y(i) << " " << U(i) << endl;
  }
  flux.close();
}

int Mesh_Adapt::cellule(double distance) //A changer surement avec l'adaptation de maillage !!!!!!!!!!!!!!!!!!
{
  int Nx=_df->Get_Nx(), Ny=_df->Get_Ny();
  double dy = _df->Get_dy();
  int Nydist = (Nx*(Ny-int(ceil(distance/dy)))-1);

  return Nydist;
}

#define _MESHADAPT_CPP
#endif
