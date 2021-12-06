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
  Eigen::VectorXd _sol;
  // Time
  double _t;

public:
  // Constructeur par défaut
  TimeScheme(DataFile* data_file, FiniteVolume* fin_vol);
  // Destructeur par défaut - Si la classe ne contient pas de destructeur par défaut
  // alors le compilateur en génère un implicitement.
  virtual ~TimeScheme();
  void InitialCondition();
  // Enregistre la solution un fichier
  void SaveSol(Eigen::VectorXd sol, std::string n_sol, int n);
  // Enregistre la solution un fichier
  // void Save_sol(Eigen::VectorXd s, int n, std::string st) {_fin_vol->Save_sol(s, n, st);};??
  // Une étape du schéma en temps
  virtual void Advance() = 0;
  // Permet de récupérer _sol
  const Eigen::VectorXd & GetSolution() const;
};

class EulerScheme : public TimeScheme
{
public:
  EulerScheme(DataFile* data_file, FiniteVolume* fin_vol);
  void Advance();
};

class ImplicitEulerScheme : public TimeScheme
{
private:
	Eigen::SparseMatrix<double> _IdplusdtsigmaH;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > _solver_direct;
public:
  ImplicitEulerScheme(DataFile* data_file, FiniteVolume* fin_vol);
  void Advance();
};

#define _TIME_SCHEME_H
#endif
