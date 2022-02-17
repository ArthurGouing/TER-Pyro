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
  Eigen::VectorXd _T;
  Eigen::VectorXd _rho;
  Eigen::VectorXd _rhostar;
 public:
   //Constructeur
   Solution();
   //Valeur aux mailles des grandeurs (surcharge de fonction)
   double T(int i, int j);
   double T_case(double x, double y);
   double T(double n);

   //sous forme de vecteur
   Eigen::VectorXd Get_T(){return _T;};  //le vecteur de taille Nx*Ny
   Eigen::VectorXd Get_Ty(int i); //le vecteur de taille Ny ...
   Eigen::VectorXd Get_Tx(int j); //le vecteur de taille Nx ...

   //idem avec rho
   double rho(int i, int j);
   double rho_case(double x, double y);
   double rho(double n);

   Eigen::VectorXd Get_rho(){return _rho;};
   Eigen::VectorXd Get_rhoy(int i); //le vecteur de taille Ny ...
   Eigen::VectorXd Get_rhox(int j); //le vecteur de taille Nx ...


   //et rho* ?
   double rhostar(int i, int j);
   double rhostar_case(double x, double y);
   double rhostar(double n);

   Eigen::VectorXd Get_rhostar(){return _rhostar;};
   Eigen::VectorXd Get_rhostary(int i); //le vecteur de taille Ny ...
   Eigen::VectorXd Get_rhostarx(int j); //le vecteur de taille Nx ...
};

#define _SOLUTION_H
#endif
