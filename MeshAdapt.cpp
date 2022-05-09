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
    Maillage_non_uniforme();
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
  double metmax=10;
  double minimum =0.01;


  vitesse();
  //U2=Derive_y_2(_rho)/5000; // dérivée seconde selon y en x = dx, aux noeuds du maillaage
  U2=Derive_y_2(sol.Get_rho()/500); // dérivée seconde selon y en x = dx, aux noeuds du maillaage

  for (int i=0; i<U2.size() ; i++)
  {
    U2(i)=exp(-pow(_Y(i)-0.005,2)/2e-7);
  }
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
    metric(i)= sqrt(max(abs(U2(i)),minimum));
    //cout << i<<" U2  ; " << U2(i) << endl;
    if (_df->Get_scenario() == "noadapt")
    {metric(i) = 1;} // A enlevé pour acitvé l'adapttion
    cout<< "metrice" << metric(i) << endl;
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
  double coeff = 1e12;
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
  fstream file;
  string a="mesh.dat";
  file.open(a, ios::out);
  cout << "file open : "<<file.is_open()<<endl;
  file<< "i ; _Y(i) ; _Yold ; _Dy ; _Dyold"<<endl;
  for (int i=1; i<_Y.size()-1;i++) //On remplit _Y
  {
    _Y(i)=Ysolv(i-1);
    file << i << " "<<_Y(i) << " " << _Yold((i))<<" "<<_Dy(i)<< " "<< _Dyold(i)<< endl;
  }
  file.close();

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
  Eigen::VectorXd U2;
  U2=Derive_y_2(sol.Get_rho()); // dérivée seconde selon y en x = dx, aux noeuds du maillaage
  int Ny = _df->Get_Ny();
  fstream file;
  fstream file1;
  string a="Mesh/rho_"+text+".dat";
  string b="Mesh/derive2_"+text+".dat";
  file.open(a, ios::out);

  file1.open(b, ios::out);

  for (int i=0; i<U2.size() ; i++)
  {
    U2(i)=exp(-pow(_Y(i)-0.005,2)/2e-7);
  }

  for(int i=0; i<_Y.size() ;i++)
  {
  //   cout <<_Y(i)<< endl;
  //   cout << sol.Get_rhoy(0)(i) << e0ndl;
  U2(i)=exp(-pow(_Y(i)-0.005,2)/2e-7);
  file1 << _Y(i) << " " << U2(i) << endl;
  U2(i)=exp(-pow(_Yold(i)-0.005,2)/2e-7);
  file << _Yold(i) <<" "<< U2(i) << endl;
}
file1 << _Y(_Y.size()-1) << " " << U2(_Y.size()-1) << endl;

for (int i=0; i<Ny+1 ; i++)
{
  // cout << endl<<"itération: "<<i<<" sur "<<Ny<<endl;
  int i_rho=-1;
  for (int k=0 ; k<_Dy.size()-1 ; k++)
  {

    if ( (_Y(i)>_Yold(k)+_Dyold(k)/2)&&(_Y(i)<_Yold(k+1)+_Dyold(k+1)/2) )
    {
      // cout <<_Yold(k+1)<< " et "<< _Dyold(k+1)<<endl;
      // cout <<_Yold(k)+_Dyold(k)/2<<"<"<<_Y(i)<<"<"<<_Yold(k+1)+_Dyold(k+1)/2<<endl;
      i_rho=k;
      // cout << i_rho<<endl;
    }
    //cout << k <<" "<<"taille de Dy : " <<_Dy.size()-2<<endl;
  }
  if (_Y(i)<_Yold(0)+_Dyold(0)/2)
  {
    file << _Y(i) << " " << sol.Get_rhoy(0)(0)<<endl;
    i_rho = i_rho-1;
    // cout << "test" << endl;
  }
  else if (_Y(i)>_Yold(_Y.size()-2)+_Dyold(_Y.size()-2)/2)
  {
    file << _Y(i) << " " << sol.Get_rhoy(0)(_Y.size()-2)<<endl;
    i_rho = i_rho-1;
    // cout << "test1" << endl;

  }
  else if (i_rho>=0) // != de 0 ou -1
  {
    double a=_Yold(i_rho)+_Dyold(i_rho)/2;
    double b=_Yold(i_rho+1)+_Dyold(i_rho+1)/2;
    double c=_Y(i);
    // cout << "On a i_rho : "<<i_rho<<endl;
    // cout << "rho(i_rho) : " << sol.Get_rhoy(0)(i_rho)<< " "<< "rho(i_rho+1) : " << sol.Get_rhoy(0)(i_rho+1)<<endl;
    // cout << "le rho sur _Y(i) : "<< ((b-c)*sol.Get_rhoy(0)(i_rho) + (c-a)*sol.Get_rhoy(0)(i_rho+1))/(b-a) <<endl;
    // cout << "le dx de la formule : "<< (c-a) <<endl;
     file<< _Y(i) << " "<< ((b-c)*sol.Get_rhoy(0)(i_rho) + (c-a)*sol.Get_rhoy(0)(i_rho+1))/(b-a) <<endl;//
  }
  else
  {
    cout <<i_rho<<endl;
    cout <<"Il y a une couille quelques parts !"<<endl;
  }
}
file.close();
file1.close();
}


//--------------------------------------Dérivée seconde-------------------------

VectorXd Mesh_Adapt:: Derive_y_2(VectorXd T) // T de taille Nx*Ny
{

  int Ny = _df->Get_Ny();
  int Nx = _df->Get_Nx();
  VectorXd derive2 (_Dy.size()+1);

  for (int i=0; i<_Dy.size();i++)
  {
    _Dy(i)=_df->Get_dy();
  }

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
void Mesh_Adapt::Maillage_non_uniforme() // Calcul du maillage non uniforme !!!!!!!
{
  double a , Ny=_df->Get_Ny() , dy=_df->Get_dy() , taille=Ny*dy , length(0);
  a=(taille-dy/2.*Ny)/(Ny*Ny-Ny*(Ny+1)/2.);
  for (int i=0 ; i<Ny ;  i++ )
  {
    _Dy(i)=(Ny-i-1)*a+dy/2.;
    cout << _Dy(i) << endl;
    length+=_Dy(i);
  }
  cout <<"la longueur totale suivant y est " << length << endl;
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
    _v(i)=(_Y(i+1)-_Y(i))/(2*dt)-(_Yold(i+1)-_Yold(i))/(2*dt); // vi+1/2
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
  for (int i=0 ; i<_Yold.size() ; i++)
  {
    _Yold(i)=i*_df->Get_dy();
  }
    //_Yold=_Y; /////////////////////////////////////////////////////
}

void Mesh_Adapt::Update_Dyprevious()
{
  _Dy = _Dyold;
  _Dyprevious=_Dy;
}


#define _MESHADAPT_CPP
#endif
