#ifndef _SOLUTION_CPP

#include "Solution.h"

Solution::Solution(data_file):
df(data_file)
{
  _T.resize(_df->Get_Nx()) * _df->Get_Ny())
  _rho.resize(_df->Get_Nx()) * _df->Get_Ny())
  _rhostar.resize(_df->Get_Nx()) * _df->Get_Ny())
}

Solution :: doouble T(double n)
{
  return _T(n);
}


































Solution::Eigen::VectorXd Get_Ty(int j)//Le vecteur de la colonne j
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int i=0; i<V.size(); i++)
  V(i)=T(i,j);
  return V;
}

Solution::Eigen::VectorXd Get_Tx(int i)// Le vecteur de la ligne i
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int j=0; j<V.size(); j++)
  V(j)=T(i,j);
  return V;
}



Solution::Eigen::VectorXd Get_rhoy(int j)//Le vecteur de la colonne j
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int i=0; i<V.size(); i++)
  V(i)=rho(i,j);
  return V;
}

Solution::Eigen::VectorXd Get_rhox(int i)// Le vecteur de la ligne i
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int j=0; j<V.size(); j++)
  V(j)=rho(i,j);
  return V;
}



Solution::Eigen::VectorXd Get_rhostary(int j)//Le vecteur de la colonne j
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int i=0; i<V.size(); i++)
  V(i)=rhostar(i,j);
  return V;
}

Solution::Eigen::VectorXd Get_rhostarx(int i)// Le vecteur de la ligne i
{
  Eigen::VectorXd V;
  V.resize(_df->Get_Ny());
  for (int j=0; j<V.size(); j++)
  V(j)=rhostar(i,j);
  return V;
}
#define _SOLUTION_CPP
#endif
