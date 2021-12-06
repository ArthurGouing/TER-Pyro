#ifndef _FUNCTION_H

#include "DataFile.h"

class Function {
private:
  // Quelques variables privées utiles pour
  // construire la condition initiale et la solution exacte (sigma)

public:
  // Constructeur
  Function(DataFile* data_file);

  // Condition initiale
  double InitialCondition(const double x, const double y) const;

  // Terme source
  double SourceFunction(const double t) const;
};

#define _FUNCTION_H
#endif
