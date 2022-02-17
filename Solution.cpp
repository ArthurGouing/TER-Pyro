#ifndef _SOLUTION_CPP

#include "Solution.h"

Solution::Solution(DataFile* data_file):
_df(data_file)
{
  T.resize(_df->Get_Nx() * _df->Get_Ny());
  rho.resize(_df->Get_Nx() * _df->Get_Ny());
  rhostar.resize(_df->Get_Nx() * _df->Get_Ny());
}

double Solution::T_cell(int i, int j) //i la ligne et j la colonne
{
  double Nx=_df->Get_Nx();
  return T[i*Nx+j];
}


double Solution::T_cell(double x, double y)
{
  double dx=_df->Get_dx();
  double dy=_df->Get_dy();
  int i=int(x/dx)+1;
  int j=int(y/dy)+1;
  return T_cell(i,j);
}

double Solution::rho_cell(int i, int j) //i la ligne et j la colonne
{
  double Nx=_df->Get_Nx();
  return rho[i*Nx+j];
}

double Solution::rho_cell(double x, double y)
{
  double dx=_df->Get_dx();
  double dy=_df->Get_dy();
  int i=int(x/dx)+1;
  int j=int(y/dy)+1;
  return rho_cell(i,j);
}


double Solution::rhostar_cell(int i, int j) //i la ligne et j la colonne
{
  double Nx=_df->Get_Nx();
  return rhostar[i*Nx+j];
}

double Solution::rhostar_cell(double x, double y)
{
  double dx=_df->Get_dx();
  double dy=_df->Get_dy();
  int i=int(x/dx)+1;
  int j=int(y/dy)+1;
  return rhostar_cell(i,j);
}


Eigen::VectorXd Solution::Get_Ty(int j)//Le vecteur de la colonne j
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int i=0; i<V.size(); i++)
  V(i)=T_cell(i,j);
  return V;
}

Eigen::VectorXd Solution::Get_Tx(int i)// Le vecteur de la ligne i
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int j=0; j<V.size(); j++)
  V(j)=T_cell(i,j);
  return V;
}



Eigen::VectorXd Solution::Get_rhoy(int j)//Le vecteur de la colonne j
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int i=0; i<V.size(); i++)
  V(i)=rho_cell(i,j);
  return V;
}

Eigen::VectorXd Solution::Get_rhox(int i)// Le vecteur de la ligne i
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int j=0; j<V.size(); j++)
  V(j)=rho_cell(i,j);
  return V;
}



Eigen::VectorXd Solution::Get_rhostary(int j)//Le vecteur de la colonne j
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int i=0; i<V.size(); i++)
  V(i)=rhostar_cell(i,j);
  return V;
}

Eigen::VectorXd Solution::Get_rhostarx(int i)// Le vecteur de la ligne i
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int j=0; j<V.size(); j++)
  V(j)=rhostar_cell(i,j);
  return V;
}

#define _SOLUTION_CPP
#endif
