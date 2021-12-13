#ifndef _MESHADAPT_H

#include "DataFile.h"
#include "Dense"
#include "Sparse"



class Mesh_Adapt
 {
 private:
   Eigen::VectorXd _T;
  Eigen::VectorXd _Dy;
  Eigen::VectorXd _Y;
  DataFile* _df;
 public:
  Mesh_Adapt(DataFile* data_file);
  void Update(Eigen::VectorXd T);
  Eigen::VectorXd Derive_y_2(Eigen::VectorXd T);
  void save_vector(Eigen::VectorXd Y, std::string a);
  const Eigen::VectorXd Get_Dy() const {return _Dy;}; /////!!!!!!
};

#define _MESHADAPT_CPP
#endif
