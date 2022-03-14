#ifndef _MESHADAPT_CPP


#include "MeshAdapt.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen"

using namespace Eigen;
using namespace std;

//----------------------------Constructeur--------------------------------------

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

}

//--------------------------------Mise à jour du maillage-----------------------

void Mesh_Adapt::Update(Solution & sol)//Soluiton sol en entrée
{
  //calcul u

  int Nx = _df->Get_Nx();
  int Ny = _df->Get_Ny();
  double Ly = _df->Get_ymax();
  VectorXd U2(Ny+1);
  VectorXd metric(Ny+1);
  VectorXd K(Ny);
  _rho= sol.Get_rho(); // inutile ?
  double maxU2 =0;
  double metmax=10.0;

  U2=Derive_y_2(_rho)/5000; // dérivée seconde selon y en x = dx, aux noeuds du maillaage
  for (int i=0; i<U2.size(); i++)
  {
    if (abs(U2(i))>maxU2)
    {
      maxU2=abs(U2(i));
    }
  }
  U2=(U2/maxU2)*metmax;

  //calcul des coefficient de raideur ki
  for (int i=0 ; i<metric.size(); i++)  //Calcul de la métrique
  {
    double minimum =10;
    metric(i)= sqrt(max(abs(U2(i)),minimum));
  }

  for (int i=0 ; i<K.size(); i++)  //Calcul de K
  {
    K(i)=(metric(i)+metric(i+1))/2.;
  }


  //On construit M et b puis on résout M.x = b
  SparseMatrix<double> M(Ny-1, Ny-1);
  VectorXd b(Ny-1), Ysolv(Ny-1);
  vector<Triplet<double>> triplets;
  SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
  double coeff = 1e10;
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

  for (int i=1; i<_Y.size()-1;i++) //On remplit M
  {
    _Y(i)=Ysolv(i-1);
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
}

//--------------------------------------Dérivée seconde-------------------------

VectorXd Mesh_Adapt:: Derive_y_2(VectorXd T) // T de taille Nx*Ny
{

  int Ny = _df->Get_Ny();
  int Nx = _df->Get_Nx();
  VectorXd derive2 (_Dy.size()+1);
  for (int i=2;i<derive2.size()-2; i++) //Cas général
  {
    double Tim1 = T[(i-1)*Nx];
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

//----------------------------Numéro de la case--------------------------------

int Mesh_Adapt::cellule(double distance) //A changer surement avec l'adaptation de maillage !!!!!!!!!!!!!!!!!!
{
  int Nx=_df->Get_Nx(), Ny=_df->Get_Ny();
  double dy = _df->Get_dy();
  int Nydist = (Nx*(Ny-int(ceil(distance/dy)))-1);

  return Nydist;
}

#define _MESHADAPT_CPP
#endif
