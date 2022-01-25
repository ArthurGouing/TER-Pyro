
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
  _Y(0)=0;
  for (int i=0; i<_Dy.size(); i++)
  {
    _Dy(i) =_df->Get_dy();
    _Y(i+1) = _Y(i)+_df->Get_dy();
  }

}

void Mesh_Adapt::Update(VectorXd rho)
{
  //calcul u
  int Nx = _df->Get_Nx();
  int Ny = _df->Get_Ny();
  double Ly = _df->Get_ymax();
  VectorXd U2(Ny+1);
  VectorXd metric(Ny+1);
  VectorXd K(Ny);
  _rho= rho;

  // for (int i=0; i<Ny+1;i++)
  // {
  //   //double nu = 2;
  //   //double c = 1;
  //   //U2(i)=c*exp(c*_Y(i)/(nu-1))/((nu-1)*(1-exp(c/nu)));
  //   U2(i)=25*exp(-(_Y(i)-0.005)*(_Y(i)-0.005)/(2*0.0002*0.0002));
  // }
  // //---------------Affichage test ------------------------------------------------
  //
  //   //On trace la D2, sur un maillage uniforme
  //   save_vector(U2,_Y, "test_d2T0.dat");
  //   //On trace la valeur de la Température
  //   save_vector(T,_Y, "test_T.dat");// ne marche pas car T = Ty
  //   //cout << "Dy : " << _Dy<< endl;
  //
  // //------------------------------------------------------------------------------
  U2=Derive_y_2(_rho); // dérivée seconde selon y en x = dx, aux noeuds du maillaage



  //calcul des coefficient de raideur ki
  for (int i=0 ; i<metric.size(); i++)  //Calcul de la métrique
  {
    double minimum = 0.01;
    metric(i)= sqrt( max(U2(i),minimum));
  }

  for (int i=0 ; i<K.size(); i++)  //Calcul de K
  {
    K(i)=(metric(i+1)+metric(i))/2;
  }

  //On construit M et b puis on résout M.x = b
  SparseMatrix<double> M(Ny-1, Ny-1);
  VectorXd b(Ny-1), Ysolv(Ny-1);
  vector<Triplet<double>> triplets;
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
  double coeff = 1;
  b.setZero();
  b(b.size()-1)=coeff*K(K.size()-1)*Ly;

  //On remplit M
  //i = 1
  triplets.push_back({0,0,coeff*(K(1)+K(0)) });
  triplets.push_back({0,1,-coeff*K(1)});

  //pour i=2 à Ny-2
  for (int i=1; i<M.rows()-1;i++)
  {
    int ischema = i+1;
    triplets.push_back({i,i-1, -K(ischema-1)});
    triplets.push_back({i,i,  +(K(ischema-1)+K(ischema)) });
    triplets.push_back({i,i+1, -K(ischema)});
  }

  //pour i=Ny-1  (M.row() = Ny-1)
  triplets.push_back({M.rows()-1,M.rows()-2,  -coeff*K(Ny-2)});
  triplets.push_back({M.rows()-1,M.rows()-1,  +coeff*(K(Ny-2)+K(Ny-1)) });

  M.setFromTriplets(triplets.begin(), triplets.end());

  solver.analyzePattern(M);
  solver.factorize(M);
  Ysolv = solver.solve(b);

  //on actualise les _Y à l'intérieur du maillage
  for (int i=0; i<Ysolv.size();i++)
  {
    _Y(i+1)= Ysolv(i);
  }


  // On détermine les _Dy a partir des nouveaux _Y
  for (int i=0;i<_Dy.size();i++)
  {
    _Dy(i)= _Y(i+1)-_Y(i);
    if (_Dy(i)==0)
    {
      cout<< "Le dy de la maille "<<i<<" vaut 0, on aura une division par 0 dans la construction des matrices"<<endl;

    }
  }

  //
  // // test de la fonction derive2 pour un maillage non uniforme
  // cout << "debut test derive" << endl;
  // VectorXd F(Ny+1), F2(Ny+1), Fexact(Ny+1);
  // for (int i=0 ; i<Ny+1 ; i++)
  // {
  //   F(i) = pow(_Y(i),3);
  //   Fexact(i) = 6*_Y(i);
  // }
  // F2 = Derive_y_2(F);
  // cout << "save"<<endl;
  // save_vector(F, _Y, "fonction.dat");
  // save_vector(Fexact, _Y, "fonctionE.dat");
  // save_vector(F2, _Y, "fonction2.dat");
  //
}

VectorXd Mesh_Adapt:: Derive_y_2(VectorXd T)
{

  int Ny = _df->Get_Ny();
  int Nx = _df->Get_Nx();
  VectorXd derive2 (_Dy.size()+1);
  for (int i=2;i<derive2.size()-2; i++) //Cas général
  {

    // double Tim1 = (T[(i-1)*Nx]+T[(i)*Nx])/2.;          // On fait 3 calculs on pourrait en faire qu'1 et récupérer les valeurs déjà calculé
    // double Ti = (T[(i)*Nx]+T[(i+1)*Nx])/2.;
    // double Tip1 = (T[(i+2)*Nx]+T[(i+1)*Nx])/2.;

    double Tim1 = T[(i-1)*Nx];                            // On fait 3 calculs on pourrait en faire qu'1 et récupérer les valeurs déjà calculé
    double Ti =  T[i*Nx];
    double Tip1 =  T[(i+1)*Nx];


    derive2(i) = 2*(_Dy(i-1)*Tip1 + _Dy(i)*Tim1 - (_Dy(i-1)+_Dy(i))*Ti) / (_Dy(i)*_Dy(i)*_Dy(i-1) + _Dy(i-1)*_Dy(i-1)*_Dy(i));
  }

  //Conditions aux bords
  derive2(0) = derive2(2);
  derive2(1) = derive2(2);
  derive2(derive2.size()-1) = derive2(derive2.size()-3);
  derive2(derive2.size()-2) = derive2(derive2.size()-3);
  return derive2;
}



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
  cout << "Dans le fichier "+a+" le vecteur U est de taill "<<U.size()<<" et Y est de taille "<<Y.size()<<endl;
  ofstream flux;
  flux.open(a);
  for (int i=0; i<Y.size();i++)
  {
    flux << Y(i) << " " << U(i) << endl;
  }

  flux.close();
}

#define _MESHADAPT_CPP
#endif
