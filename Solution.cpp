#ifndef _SOLUTION_CPP

#include "Solution.h"

Solution::Solution(data_file):
df(data_file)
{
<<<<<<< HEAD
  _T.resize(df->Get_Nx()) * df->Get_Ny());
  _rho.resize(df->Get_Nx()) * df->Get_Ny());
  _rhostar.resize(df->Get_Nx()) * df->Get_Ny());
=======
  _T.resize(_df->Get_Nx()) * _df->Get_Ny())
  _rho.resize(_df->Get_Nx()) * _df->Get_Ny())
  _rhostar.resize(_df->Get_Nx()) * _df->Get_Ny())
>>>>>>> e03658955170f057c746233fe8a191fb94400ae1
}

Solution :: doouble T(double n)
{
  return _T(n);
<<<<<<< HEAD
}

Solution :: double T(int i, int j) //i la ligne et j la colonne
{
  Nx=_df->Get_Nx();
  return _T(i*Nx+j);
}


Solution :: double T_case(double x, double y)
{
  dx=_df->Get_dx();
  dy=_df->Get_dy();
  i=int(x/dx)+1;
  j=int(y/dy)+1;
  return T(i,j);
}


Solution :: double rho(double n)
{
  return _rho(n);
}


Solution :: double rho(int i, int j) //i la ligne et j la colonne
{
  Nx=_df->Get_Nx();
  return _rho(i*Nx+j);
}

Solution :: double rho_case(double x, double y)
{
  dx=_df->Get_dx();
  dy=_df->Get_dy();
  i=int(x/dx)+1;
  j=int(y/dy)+1;
  return rho(i,j);
=======
>>>>>>> e03658955170f057c746233fe8a191fb94400ae1
}


Solution :: double rhostar(double n)
{
  return _rhostar(n);
}


Solution :: double rhostar(int i, int j) //i la ligne et j la colonne
{
  Nx=_df->Get_Nx();
  return _rhostar(i*Nx+j);
}

Solution :: double rhostar_case(double x, double y)
{
  dx=_df->Get_dx();
  dy=_df->Get_dy();
  i=int(x/dx)+1;
  j=int(y/dy)+1;
  return rhostar(i,j);
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
