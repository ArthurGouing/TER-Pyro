#ifndef _SOLUTION_H

/*la class Solution est utilisé dans les fichier functions ; ou pas
J'ai pas modifié les _sol et _rho dans la partie Euler Implicit de TimeScheme
on verra comment on fait avec _rhostar
le reste devrait être bon*/

#include "DataFile.h"
#include "Dense"
#include "Sparse"


class Solution
 {
 private:
  Eigen::VectorXd _T;
  Eigen::VectorXd _rho;
 public:
   //Constructeur
   Solution();
   //Valeur aux mailles des grandeurs (surcharge de fonction)
   double T(int i, int j);
   double T(double x, double y);
   double T(double n);

   //sous forme de vecteur
   Eigen::VectorXd Get_T(){return _T;};  //le vecteur de taille Nx*Ny
   Eigen::VectorXd Get_Ty(); //le vecteur de taille Ny ...
   Eigen::VectorXd Get_Tx(); //le vecteur de taille Nx ...

   //idem avec rho
   //et rho* ?






};

#define _SOLUTION_H
#endif
