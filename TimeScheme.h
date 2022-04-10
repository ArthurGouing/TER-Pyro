#ifndef _TIME_SCHEME_H

#include "FiniteVolume.h"

class TimeScheme
{
protected:
  // Pointeur vers le fichier de données
  DataFile* _df;
  // Pointeur vers le laplacien
  FiniteVolume* _fin_vol;
  //Vecteur initial et vecteur solution
  //Eigen::VectorXd _sol, _rho, _rhostar;
  Solution _sol, _solold; //_soldold contient ce dont on a besoin comme CI pour faire l'adaptation de maillage
  // Time
  double _t;

public:
  // Constructeur par défaut
  TimeScheme(DataFile* data_file, FiniteVolume* fin_vol);
  // Destructeur par défaut - Si la classe ne contient pas de destructeur par défaut
  // alors le compilateur en génère un implicitement.
  virtual ~TimeScheme();
  void InitialCondition();
  void Update_Told_rhoold();
  // Enregistre la solution un fichier
  void SaveSol(Solution sol, std::string n_sol, int n);
  //void SaveSol(Eigen::VectorXd sol, std::string n_sol, int n); //à changer
  void Save_rho(Eigen::VectorXd rho , double t , std::string name_file);//à changer
  // Une étape du schéma en temps
  virtual void Advance() = 0;
  virtual void Advance_ALE() = 0;
  // Permet de récupérer _sol
  /*
  const Eigen::VectorXd & GetSolution() const;
  const Eigen::VectorXd & GetSolutionrho() const; // Inutile à cause de la class S
  */
  Solution & Get_Solution(){return _sol;};
};


class ImplicitEulerScheme : public TimeScheme
{
private:
  Eigen::SparseLU<Eigen::SparseMatrix<double> > _solver_direct;
public:
  ImplicitEulerScheme(DataFile* data_file, FiniteVolume* fin_vol);
  void Advance();
  void Advance_ALE();
};

#define _TIME_SCHEME_H
#endif
