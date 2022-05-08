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
  _Dystar.resize(_df->Get_Ny());
  _Dyold.resize(_df->Get_Ny());
  _Dyprevious.resize(_df->Get_Ny());
  _Y.resize(_df->Get_Ny()+1);
  _Yold.resize(_df->Get_Ny()+1);
  _v.resize(_df->Get_Ny()+1);
  _Y(0)=0.;
  if (_df->Get_CastestnonUnif() == "oui") //ajout du cas maillage non uniforme  attention au datafile!!!!!!!
  {
    _Dy=Maillage_non_uniforme();
    _Dyold=Maillage_non_uniforme();
  }
  else
  {
    for (int i=0; i<_Dy.size(); i++)
    {
      _Dy(i) =_df->Get_dy();
      _Dyold(i)=_df->Get_dy();
      _Dyprevious(i) = _df->Get_dy();
      _Y(i+1) = _Y(i)+_df->Get_dy();
    }
  }
}

//--------------------------------Mise à jour du maillage-----------------------

void Mesh_Adapt::Update(Solution & sol)
//  Enlevé la ligne 76 pour réactiver l'adaptation de maillage
{

  //Initialisation
  int Nx = _df->Get_Nx();
  int Ny = _df->Get_Ny();
  double Ly = _df->Get_ymax();
  VectorXd U2(Ny+1);
  VectorXd metric(Ny+1);
  VectorXd K(Ny);
  //_rho= sol.Get_rho(); // inutile ?
  double maxU2 =1;
  double metmax=15;

  vitesse();
  //U2=Derive_y_2(_rho)/5000; // dérivée seconde selon y en x = dx, aux noeuds du maillaage
  U2=Derive_y_2(sol.Get_rho())/5000; // dérivée seconde selon y en x = dx, aux noeuds du maillaage
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
    double minimum =1;
    metric(i)= sqrt(max(abs(U2(i)),minimum));
    //metric(i) = 1; // A enlevé pour acitvé l'adapttion
  }

  for (int i=0 ; i<K.size(); i++)  //Calcul de K
  {
    K(i)=(metric(i)+metric(i+1))/2.;
    if (K(i)==0)
    {
      cout<<"K=0"<< endl;
    }
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

cout <<"4"<<endl;
  Ysolv = solver.solve(b);

  cout <<"5"<<endl;

//Pour le test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  _Dy=Maillage_non_uniforme();
  fstream file;
  string a="mesh.dat";
  file.open(a, ios::out);
  for (int i=1; i<_Y.size()-1;i++) //On remplit _Y
  {
    _Y(i)=Ysolv(i-1); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!vrai truc
    file << _Y(i) << " " << U2(i) << endl;
  }
  file.close();
  ///////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // On détermine les _Dy a partir des nouveaux _Y calculés
  for (int i=0;i<_Dy.size();i++)
  {
    _Dy(i)= _Y(i+1)-_Y(i);
    if (_Dy(i)<=1e-16)
    {
      cout<< "Le dy de la maille "<<i<<" vaut tres petit, on aura une division par 0 dans la construction des matrices"<<endl;

    }
  }
}

void Mesh_Adapt::Update2(Solution & sol)
//  Enlevé la ligne 76 pour réactiver l'adaptation de maillage
{

  //Initialisation
  int Nx = _df->Get_Nx();
  int Ny = _df->Get_Ny();
  double Ly = _df->Get_ymax();
  VectorXd U2(Ny+1);
  VectorXd metric(Ny+1);
  VectorXd K(Ny);
  //_rho= sol.Get_rho(); // inutile ?
  double maxU2 =1;
  double metmax=15;

//Pour le test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  _Dy=Maillage_non_uniforme();
  fstream file;
  string a="mesh.dat";
  file.open(a, ios::out);
  for (int i=1; i<_Y.size()-1;i++) //On remplit _Y
  {
    //_Y(i)=Ysolv(i-1); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!vrai truc
    _Y(i)=_Y(i-1)+_Dy(i-1); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    file << _Y(i) << " " << U2(i) << endl;
  }
  file.close();
  ///////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // On détermine les _Dy a partir des nouveaux _Y calculés
  for (int i=0;i<_Dy.size();i++)
  {
    _Dy(i)= _Y(i+1)-_Y(i);
    if (_Dy(i)<=1e-16)
    {
      cout<< "Le dy de la maille "<<i<<" vaut tres petit, on aura une division par 0 dans la construction des matrices"<<endl;

    }
  }
}

void Mesh_Adapt::Update3(Solution & sol)
//  Enlevé la ligne 76 pour réactiver l'adaptation de maillage
{

  //Initialisation
  int Nx = _df->Get_Nx();
  int Ny = _df->Get_Ny();
  double Ly = _df->Get_ymax();
  VectorXd U2(Ny+1);
  VectorXd metric(Ny+1);
  VectorXd K(Ny);
  //_rho= sol.Get_rho(); // inutile ?
  double maxU2 =1;
  double metmax=15;

//Pour le test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  _Dy=Maillage_non_uniforme2();
  fstream file;
  string a="mesh.dat";
  file.open(a, ios::out);
  for (int i=1; i<_Y.size()-1;i++) //On remplit _Y
  {
    //_Y(i)=Ysolv(i-1); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!vrai truc
    _Y(i)=_Y(i-1)+_Dy(i-1); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    file << _Y(i) << " " << U2(i) << endl;
  }
  file.close();
  ///////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // On détermine les _Dy a partir des nouveaux _Y calculés
  for (int i=0;i<_Dy.size();i++)
  {
    _Dy(i)= _Y(i+1)-_Y(i);
    if (_Dy(i)<=1e-16)
    {
      cout<< "Le dy de la maille "<<i<<" vaut tres petit, on aura une division par 0 dans la construction des matrices"<<endl;

    }
  }
}

void Mesh_Adapt::Affichage(std::string text,Solution sol)
{
  fstream file;
  string a="Mesh/mesh"+text+".dat";
  cout <<a << endl;
  file.open(a, ios::out);
  for(int i=0; i<_Y.size() ;i++)
  {
    file << _Y(i) << " "<< 0 << endl;
  }
  file.close();
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


int Mesh_Adapt::cellule(double distance)
{
  int Nx=_df->Get_Nx(), Ny=_df->Get_Ny();
  int Nydist; //(Nx*(Ny-int(ceil(distance/dy)))-1);
  double sommedist = 0., sommedist2=0.;
  int i=0;

  while (sommedist<=distance)
  {
    sommedist+=_Dy(Ny-1-i);
    i+=1;
  }
  sommedist2=sommedist-_Dy(Ny-1-(i-1));
  if (abs(sommedist2-distance)<abs(sommedist-distance))
  {
    Nydist = (Nx*(Ny-(i-1))-1);
  }
  else
  {
    Nydist = (Nx*(Ny-i)-1);
  }

  return Nydist;
}

//----------------------------------------------------------------------------
Eigen::VectorXd Mesh_Adapt::Maillage_non_uniforme() // Calcul du maillage non uniforme !!!!!!!
{
  VectorXd Dy;
  Dy.resize(_df->Get_Ny());
  double a , Ny=_df->Get_Ny() , dy=_df->Get_dy() , taille=Ny*dy , length(0);
  a=(taille-dy/2.*Ny)/(Ny*Ny-Ny*(Ny+1)/2.);
  for (int i=0 ; i<Ny ;  i++ )
  {
    Dy(i)=(Ny-i-1)*a+dy/2.;
    cout << Dy(i) << endl;
    length+=Dy(i);
  }
  cout <<"la longueur totale suivant y est " << length << endl;

  return Dy;
}


Eigen::VectorXd Mesh_Adapt::Maillage_non_uniforme2() // Calcul du maillage non uniforme !!!!!!!
{
  VectorXd Dy;
  Dy.resize(_df->Get_Ny());
  double a , Ny=_df->Get_Ny() , dy=_df->Get_dy() , taille=Ny*dy , length(0);
  a=(taille-dy/2.*Ny)/(Ny*Ny-Ny*(Ny+1)/2.);
  for (int i=0 ; i<Ny ;  i++ )
  {
    Dy(Ny-1-i)=(Ny-i-1)*a+dy/2.;
    cout << Dy(Ny-1-i) << endl;
    length+=Dy(Ny-1-i);
  }
  cout <<"la longueur totale suivant y est " << length << endl;

  return Dy;
}

//---------------------------------------------------------------------------
void Mesh_Adapt::Maillage_Dystar()
{
  _Dystar = (_Dy+_Dyold)/2. ;
}

//----------------------------------------------------------------------
double Mesh_Adapt::NormLinf()
{
  double norm(0.);
  for (int i=0; i < _Dy.size(); i++)
  {
    if (abs(_Dy(i)-_Dyprevious(i))>norm)
    {
      norm=abs(_Dy(i)-_Dyprevious(i));
    }
  }
  return norm;
}
//------------------------------------------------------------------------------
void Mesh_Adapt::vitesse()  ///Y(i+1/2) est le milieu entre 2 noeuds
{
  double dt = _df->Get_dt();
  _v(0)=0.;
  for(int i=1; i<_v.size()-1;i++)
  {
    _v(i)=(_Y(i)-_Yold(i))/dt; // vi+1/2
    // cout << "_Y("<<i<<")="<< _Y(i) << endl;
    // cout << "_Yold("<<i<<")="<< _Yold(i) << endl;
    // cout << "_v("<<i<<")="<<_v(i) << endl;
  }
  _v(_v.size()-1)=0.;
}
//------------------------------------------------------------------------------

void Mesh_Adapt::Update_Dystar_vitesse()
{
  vitesse();
  Maillage_Dystar();
}

void Mesh_Adapt::Update_Dyold()
{
  _Dyold=_Dy;
  /////Debug
  cout << "Update de _Dyold"<< endl;
  for(int i=0; i<_Dyold.size(); i++)
  {
    cout << "_Dyold(" << i << ")=" << _Dyold(i) << endl;
  }
  _Yold=_Y; /////////////////////////////////////////////////////
}

void Mesh_Adapt::Update_Dyprevious()
{
  _Dyprevious=_Dy;
}


#define _MESHADAPT_CPP
#endif
