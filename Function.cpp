#ifndef _FUNCTION_CPP

#include "Function.h"
#include <cmath>

#include "Dense"
#include "Sparse"

using namespace Eigen;

Function::Function(DataFile* data_file):
_df(data_file)
{
}

// Condition initiale température
double Function::InitialCondition(const double x, const double y) const
{
  double T0=_df->Get_T0();

  return T0;
}

// Condition initiale masse volumique
double Function::InitialConditionrho(const double x, const double y) const
{
  double rhov=_df->Get_rhov();

  return rhov; //rho_vierge
}

// Terme source
double Function::SourceFunction(const double t) const
{
  if (t<=50.)
  {
    return 10000.*t;
  }
  else
  {
    return 500000.-9000.*(t-50.);
  }
}

VectorXd Function::Arrhenius(VectorXd rho, VectorXd T)
{
  double Aref=_df->Get_Aref(), Ta=_df->Get_Ta(), rhov=_df->Get_rhov(), rhop=_df->Get_rhop();
  VectorXd A(T.size());
  for (int i=0; i<T.size(); ++i)
  {
    A(i)=-Aref*exp(-Ta/T(i))*rhov*((rho(i)-rhop)/(rhov-rhop));//N ordre de la méthode
  }
  return A;
}

#define _FUNCTION_CPP
#endif
