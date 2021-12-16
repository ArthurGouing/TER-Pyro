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
  int Nx = _df->Get_Nx();
  int Ny = _df->Get_Ny();
  double Ly = _df->Get_ymax(); // = Ny*_df->Get_dy();
  VectorXd U2(Ny+1);
  VectorXd metric(Ny+1); // pour pas faire 2 fois le meme calcul de metric
  VectorXd K(Ny);
  //_Dy de taille Ny
  //_Y de taill  Ny+1
  _T= T;

for (int i=0; i<Ny+1;i++)
{
  //double nu = 2;
  //double c = 1;
  //U2(i)=c*exp(c*_Y(i)/(nu-1))/((nu-1)*(1-exp(c/nu)));
  U2(i)=25*exp(-(_Y(i)-0.005)*(_Y(i)-0.005)/(2*0.0002*0.0002));
}
  //U2=Derive_y_2(_T); // dérivée seconde selon y en x = dx, aux noeud du maillaage // Attention aux bord


  //calcul de k
  for (int i=0 ; i<metric.size(); i++)  //Calcul de la métrique
  {
    double minimum = 0.1;
    metric(i)= sqrt( max(U2(i),minimum));
  }

  for (int i=0 ; i<K.size(); i++)  //Calcul de K
  {
    K(i)=(metric(i+1)+metric(i))/2;
  }

  //on construit M et b puis on résout M.x = b
  SparseMatrix<double> M(Ny-1, Ny-1);
  VectorXd b(Ny-1), Ysolv(Ny-1);
  vector<Triplet<double>> triplets;
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;

  save_vector_mesh(_Y, "test_y0.dat");
  save_vector(U2,_Y, "test_U2.dat");

  b.setZero();
  b(b.size()-1)=K(K.size()-1)*Ly;

  //On remplit M
  //i = 1
  triplets.push_back({0,0,(K(1)+K(0)) });
  triplets.push_back({0,1,- K(1)});
  //pour i=2 à Ny-2
  for (int i=1; i<M.rows()-1;i++)
  {
    triplets.push_back({i,i-1, -K(i-1)});
    triplets.push_back({i,i,  +(K(i-1)+K(i)) });
    triplets.push_back({i,i+1, -K(i)});
  }
  //M.row() = Ny-1
  triplets.push_back({M.rows()-1,M.rows()-2,  -K(M.rows()-2)});
  triplets.push_back({M.rows()-1,M.rows()-1,  +K(M.rows()-2)+K(M.rows()-1) });

  M.setFromTriplets(triplets.begin(), triplets.end());

  solver.analyzePattern(M);
  solver.factorize(M);
  Ysolv = solver.solve(b);

  //on actualise les _Y à l'intérieur du maillage
  for (int i=0; i<Ysolv.size();i++)
  {
    _Y(i+1)= Ysolv(i);
  }

  // On détermine les dy a partir de des Xi
  _Dy(0)=_Y(0);
  _Dy(_Dy.size()-1)=Ly-_Y(_Y.size()-1);
  for (int i=1;i<_Dy.size()-1;i++)
  {
    _Dy(i)= _Y(i)-_Y(i-1);
  }

  save_vector_mesh(_Y, "test_y.dat");

  cout <<"marche !!!"<<endl;
}

VectorXd Mesh_Adapt:: Derive_y_2(VectorXd T)
{
  double Ny = _df->Get_Ny();
  double Nx = _df->Get_Nx();
  VectorXd derive2 (_Dy.size()+1);
  cout << "fin"<<endl;

  for (int i=2;i<derive2.size()-2; i++) //Cas général
  {
    double Tim1 = (T[(i-2)*Nx]+T[(i-1)*Nx])/2.;
    double Ti = (T[(i-1)*Nx]+T[i*Nx])/2.;
    double Tip1 = (T[(i+1)*Nx]+T[i*Nx])/2.;

    derive2(i) = 2*(_Dy(i)*Tip1 + _Dy(i+1)*Tim1 - (_Dy(i)+_Dy(i+1))*Ti) / (_Dy(i+1)*_Dy(i+1)*_Dy(i) + _Dy(i)*_Dy(i)*_Dy(i+1)); //expression à simplifier éventuellment
  }
  cout << "fin"<<endl;

  //Conditions aux bords
  derive2(0) = derive2(2);
  derive2(1) = derive2(2);
  derive2(derive2.size()-1) = derive2(derive2.size()-2);
  derive2(derive2.size()-2) = derive2(derive2.size()-2);
  cout << "fin derive"<<endl;
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
  ofstream flux;
  flux.open(a);
  cout << "debut save"<<endl;
  for (int i=0; i<Y.size();i++)
  {
    flux << Y(i) << " " << U(i) << endl;
  }
  cout << "debut save"<<endl;

  flux.close();
}

#define _MESHADAPT_CPP
#endif
