#ifndef _FUNCTION_H

#include "DataFile.h"
#include "Dense"
#include "Sparse"

class Function
{
private:
  // Quelques variables privées utiles pour
  // construire la condition initiale et la solution exacte (sigma)
  DataFile* _df;


public:
  // Constructeur
  Function(DataFile* data_file);

  // Condition initiale
  double InitialCondition(const double x, const double y) const;

  // Terme source
  double SourceFunction(const double t) const;

  //Arrhénius
  Eigen::VectorXd Arrhenius(Eigen::VectorXd rho, Eigen::VectorXd T);
};

#define _FUNCTION_H
#endif
