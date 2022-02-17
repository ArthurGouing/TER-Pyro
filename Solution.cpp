#ifndef _SOLUTION_CPP

#include "Solution.h"

Solution::Solution()
{
  _T.resize(df->Get_Nx()) * df->Get_Ny());
  _rho.resize(df->Get_Nx()) * df->Get_Ny());
  _rhostar.resize(df->Get_Nx()) * df->Get_Ny());
}

Solution :: doouble T(double n)
{
  return _T(n);
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





























Solution::Eigen::VectorXd Get_Ty(int i)
{
  Eigen::VectorXd Ty;
  Ty.resize(_df->Get_Ny());
}


#define _SOLUTION_CPP
#endif
