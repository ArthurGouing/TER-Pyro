#ifndef _FUNCTION_CPP

#include "Function.h"
#include <cmath>

Function::Function(DataFile* data_file)
{
}

// Condition initiale
double Function::InitialCondition(const double x, const double y) const
{
  return 293.;
}

// Terme source
double Function::SourceFunction(const double t) const
{
  if (t<=50.)
  {
    return 10000*t;
  }
  else
  {
    return 500000.-9000.*(t-50.);
  }
}

#define _FUNCTION_CPP
#endif
