#ifndef _MESHADAPT_H

#include "TimeScheme.h"


class Mesh_Adapt
 {
 private:
  Eigen::VectorXd _Dy;
  Eigen::VectorXd _X;
  DataFile* _df;
  TimeScheme* _TSch;
 public:
  Mesh_Adapt(DataFile* data_file,TimeScheme* Time_Scheme);
  void Update();
  Eigen::VectorXd Derive_y_2(Eigen::VectorXd T);
  void save_vector(Eigen::VectorXd X, std::string a);
};

#define _MESHADAPT_CPP
#endif
