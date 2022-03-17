#ifndef _MESHADAPT_H

//#include "Solution.h" // En fait on en a besoin qu'à partir de timescheme
#include "DataFile.h"
#include "Dense"
#include "Sparse"



class Mesh_Adapt
 {
 private:
  Eigen::VectorXd _rho; // inutile ; peut être pas car on fait une référence vers sol
  //Vecteur contenant les taille dy de chaques mailles
  Eigen::VectorXd _Dy;
  //Vecteur contenant les positions des "noeuds" de chaques mailles selon Y
  Eigen::VectorXd _Y;
  Eigen::VectorXd _Y1;

  //Pointeur vers le fichier de données
  DataFile* _df;
 public:
   //Constructeur
  Mesh_Adapt(DataFile* data_file);
  void Update(Solution & sol);
  Eigen::VectorXd Derive_y_2(Eigen::VectorXd rho);
  const Eigen::VectorXd Get_Dy() const {return _Dy;};
  //Fonction de test
  void save_vector_mesh(Eigen::VectorXd Y, std::string a);
  void save_vector(Eigen::VectorXd U, Eigen::VectorXd Y, std::string a); // pour voir U2
  int cellule (double distance);
  void save_vector1(Eigen::VectorXd Y, Eigen::VectorXd Y1,std::string a);
  Eigen::VectorXd vitesse();
  Eigen::VectorXd moy_arithmetique();
};

#define _MESHADAPT_H
#endif
