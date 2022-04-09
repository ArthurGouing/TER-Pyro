#ifndef _MESHADAPT_H

#include "Solution.h" // En fait on en a besoin qu'à partir de timescheme
#include "DataFile.h"
#include "Dense"
#include "Sparse"



class Mesh_Adapt
 {
 private:
  Eigen::VectorXd _rho, _Told; // inutile ; peut être pas car on fait une référence vers sol
  //Vecteur contenant les taille dy de chaques mailles
  Eigen::VectorXd _Dy,_Dystar,_Dyold, _Dyprevious;
  //Vecteur contenant les positions des "noeuds" de chaques mailles selon Y
  Eigen::VectorXd _Y,_Y1; // _Y1 pas bon défini sur le mauvais maillage
  Eigen::VectorXd _v;
  //Pointeur vers le fichier de données
  DataFile* _df;
 public:
   //Constructeur
  Mesh_Adapt(DataFile* data_file);
  void Maillage_non_uniforme();
  void Update(Solution & sol);
  Eigen::VectorXd Derive_y_2(Eigen::VectorXd rho);
  const Eigen::VectorXd Get_Dy() const {return _Dy;};
  //Fonction numéro de cellule
  const Eigen::VectorXd Get_Dystar() const {return _Dystar;};
  const Eigen::VectorXd Get_Dyold()  const {return _Dyold;};
  void Maillage_Dystar();
  double NormLinf();
  int cellule (double distance);
  void Update_Dyold();
  void Update_Dystar_vitesse();
  void Update_Dyprevious();
  void vitesse();
  const Eigen::VectorXd Get_vitesse() const {return _v;}
};

#define _MESHADAPT_H
#endif
