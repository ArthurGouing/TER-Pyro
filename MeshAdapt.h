#ifndef _MESHADAPT_H

#include "DataFile.h"
#include "Dense"
#include "Sparse"



class Mesh_Adapt
 {
 private:
   Eigen::VectorXd _T, _U;
  Eigen::VectorXd _Dy;
  Eigen::VectorXd _Y;
  DataFile* _df;
 public:
  Mesh_Adapt(DataFile* data_file);
  void Update(Eigen::VectorXd T);
  Eigen::VectorXd Derive_y_2(Eigen::VectorXd T);
  void save_vector_mesh(Eigen::VectorXd Y, std::string a);
  void save_vector(Eigen::VectorXd U, Eigen::VectorXd Y, std::string a); // pour voir U2
  const Eigen::VectorXd Get_Dy() const {return _Dy;}; /////!!!!!!
  int cellule(double distance);
};

#define _MESHADAPT_CPP
#endif
