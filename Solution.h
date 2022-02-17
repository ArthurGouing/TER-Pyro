#ifndef _SOLUTION_H

/*lLes fonction Get et Set sonnt à suprrimer une fois que les fonctioins auront
été vérifié*/

#include "DataFile.h"
#include "Dense"
#include "Sparse"


class Solution
{
public:
  //Variable
  Eigen::VectorXd T;
  Eigen::VectorXd rho;
  Eigen::VectorXd rhostar;
private:
  // Pointeur vers le fichier de données
  DataFile* _df;
public:
  //Constructeur
  Solution(DataFile* data_file);
  //Valeur aux mailles des grandeurs (surcharge de fonction)
  double T_cell(int i, int j);
  double T_cell(double x, double y);

  //sous forme de vecteur
  Eigen::VectorXd Get_T(){return T;};  //le vecteur de taille Nx*Ny
  Eigen::VectorXd Get_Ty(int i); //le vecteur de taille Ny ...
  Eigen::VectorXd Get_Tx(int j); //le vecteur de taille Nx ...
  void Set_T(int n, double a){T(n) = a;}

  //idem avec rho
  double rho_cell(int i, int j);
  double rho_cell(double x, double y);

  Eigen::VectorXd Get_rho(){return rho;};
  Eigen::VectorXd Get_rhoy(int i); //le vecteur de taille Ny ...
  Eigen::VectorXd Get_rhox(int j); //le vecteur de taille Nx ...
  void Set_rho(int n, double a){rho(n) = a;}


  //et rho* ?
  double rhostar_cell(int i, int j);
  double rhostar_cell(double x, double y);

  Eigen::VectorXd Get_rhostar(){return rhostar;};
  Eigen::VectorXd Get_rhostary(int i); //le vecteur de taille Ny ...
  Eigen::VectorXd Get_rhostarx(int j); //le vecteur de taille Nx ...
  void Set_rhostar(int n, double a){rhostar(n) = a;}
};

#define _SOLUTION_H
#endif
