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

void Mesh_Adapt::Update(VectorXd T)
{

  //calcul u
  double Nx = _df->Get_Nx();
  double Ny = _df->Get_Ny();
  double Ly= Ny*_df->Get_dy();
  VectorXd U2(_Dy.size()+1);
  _T= T;

for (int i=0; i<_Dy.size()+1;i++)
{
  U2(i)=0;
}
  //U2=Derive_y_2(_T); // dérivée seconde selon y en x = dx, aux noeud du maillaage // Attention aux bord


  //calcul de k
  VectorXd metric(_Dy.size()+1); // pour pas faire 2 fois le meme calcul de metric
  VectorXd K(_Dy.size());

  for (int i=0 ; i<metric.size(); i++)  //Calcul de la métrique
  {
    double minimum = 1;
    metric(i)= sqrt( max(U2(i),minimum));
  }

  for (int i=1 ; i<K.size(); i++)  //Calcul de K
  {
    K(i)=(metric(i+1)+metric(i))/2;
  }

  //on construit M et b puis on résout M.x = b
  SparseMatrix<double> M(_Y.size()-2, _Y.size()-2);
  VectorXd b(_Y.size()-2);
  vector<Triplet<double>> triplets;
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;

  save_vector(_Y, "test_y0");

  b.setZero();
  b(b.size()-1)=K(K.size()-1)*Ly;

  triplets.push_back({0,0,K(1)+K(0)});
  triplets.push_back({0,1,-K(0)});


  for (int i=1; i<M.rows()-1;i++) //On remplit M
  {
    triplets.push_back({i,i-1,-K(i)});
    triplets.push_back({i,i,K(i+1)+K(i)});
    triplets.push_back({i,i+1,-K(i)});
  }

  triplets.push_back({M.rows()-1,M.rows()-2,-K(M.rows()-1)});
  triplets.push_back({M.rows()-1,M.rows()-1,K(M.rows()-1)+K(M.rows())});

  M.setFromTriplets(triplets.begin(), triplets.end());

  solver.analyzePattern(M);
  solver.factorize(M);

  _Y = solver.solve(b);


  // On détermine les dy a partir de des Xi
  _Dy(0)=_Y(0);
  _Dy(_Dy.size()-1)=Ly-_Y(_Y.size()-1);
  for (int i=1;i<_Dy.size()-1;i++)
  {
    _Dy(i)= _Y(i)-_Y(i-1);
  }

  save_vector(_Y, "test_Y_update");
  save_vector(U2, "test_U2");

  cout <<"marche !!!"<<endl;
}

VectorXd Mesh_Adapt:: Derive_y_2(VectorXd T)
{
  double Ny = _df->Get_Ny();
  double Nx = _df->Get_Nx();
  VectorXd derive2 (_Dy.size()+1);

  for (int i=2;i<derive2.size()-2; i++) //Cas général
  {
    double Tim1 = (T[(i-2)*Nx]+T[(i-1)*Nx])/2.;
    double Ti = (T[(i-1)*Nx]+T[i*Nx])/2.;
    double Tip1 = (T[(i+1)*Nx]+T[i*Nx])/2.;

    derive2(i) = 2*(_Dy(i)*Tip1 + _Dy(i+1)*Tim1 - (_Dy(i)+_Dy(i+1))*Ti) / (_Dy(i+1)*_Dy(i+1)*_Dy(i) + _Dy(i)*_Dy(i)*_Dy(i+1)); //expression à simplifier éventuellment
  }

  //Conditions aux bords
  derive2(0) = derive2(2);
  derive2(1) = derive2(2);
  derive2(derive2.size()-1) = derive2(derive2.size()-2);
  derive2(derive2.size()-2) = derive2(derive2.size()-2);
  
  return derive2;
}



void Mesh_Adapt::save_vector(Eigen::VectorXd Y, std::string a) // pour voir U2 et x
{
  ofstream flux;
  flux.open(a);
  for (int i=0; i<Y.size();i++)
  {
    flux << i << " " << Y(i) << endl;
  }
  flux.close();
}



#define _MESHADAPT_CPP
#endif
