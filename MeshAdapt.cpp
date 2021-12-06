#ifndef _MESHADAPT_CPP


#include "MeshAdapt.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace Eigen;
using namespace std;


Mesh_Adapt::Mesh_Adapt(DataFile* data_file, TimeScheme* Time_Scheme) : 
  _df(data_file), _TSch(Time_Scheme)
{
}

  void Mesh_Adapt::Update()
  {
//calcul u"
double Nx = _df->Get_Nx();
VectorXd T= _TSch->GetSolution();
VectorXd U2(_Dy.size()+1);

U2=Derive_y_2(T); // dérivée seconde selon y en x = dx, aux noeud du maillaage // Attention aux bord
//calcul de k
VectorXd metric(_Dy.size()+1); // pour pas faore 2 fois le meme calcul de metric
VectorXd K(_Dy.size());
for (int i=0 ; i<metric.size(); i++)
  {
double minimum = 1;
metric(i)= sqrt( max(U2(i),minimum));
}
for (int i=0 ; i<K.size(); i++)
  {
K(i)=(metric(i-1)+metric(i))/2;
}
//on construit M et b puis on résout M.x = b
SparseMatrix<double> M(_Dy.size()-1, _Dy.size()-1);
vector<Triplet<double>> triplets;
VectorXd b(_Dy.size()-1);

b(_Dy.size()-1)=K(M.rows());
triplets.push_back({0,1,-K(0)});
triplets.push_back({0,0,K(1)+K(0)});
for (int i=1; i<M.rows()-1;i++)
  {
triplets.push_back({i,i-1,-K(i)});
triplets.push_back({i,i,K(i+1)+K(i)});
triplets.push_back({i,i+1,-K(i)});
}
triplets.push_back({M.rows(),M.rows()-1,-K(M.rows())});
triplets.push_back({M.rows(),M.rows(),K(M.rows())+K(M.rows()-1)});

M.setFromTriplets(triplets.begin(), triplets.end());

SparseLU<SparseMatrix<double>, COLAMDOrdering<int> > solver;
solver.analyzePattern(M);
solver.factorize(M);
_X = solver.solve(b);
save_vector(_X, "test_");
// On détermine les dy a partir de des Xi
_Dy(0)=_X(0);
_Dy(_Dy.size()-1)=1-_X(_X.size()-1);
for (int i=1;i<_Dy.size()-1;i++)
  {
_Dy(i)= _X(i)-_X(i-1);
}
save_vector(_X, "test_x");
save_vector(U2, "test_U2");
}

    VectorXd Mesh_Adapt:: Derive_y_2(VectorXd T)
    {

return T;
}

      void Mesh_Adapt::save_vector(Eigen::VectorXd X, std::string a) // pour voir U2 et x
      {
ofstream flux;
flux.open(a);
for (int i=0; i<X.size();i++)
  {
flux << i << " " << X(i) << endl;
}
flux.close();
}



#define _MESHADAPT_CPP
#endif
