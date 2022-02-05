#ifndef _SOLUTION_H

#include "DataFile.h"
#include "Dense"
#include "Sparse"


class Solution
 {
 private:
  Eigen::VectorXd _T;
  Eigen::VectorXd _rho;
 public:
   //Avec des surcharges de fonction
   double T(int i, int j);
   double T(double x, double y);
   double T(double n);

   //sous forme de vecteur
   Eigen::VectorXd Get_T();  //le vecteur de taille Nx*Ny
   Eigen::VectorXd Get_Ty(); //le vecteur de taille Ny
   Eigen::VectorXd Get_Tx(); //le vecteur de taille Nx

   //idem avec rho et rho* ?






};

#define _SOLUTION_H
#endif
