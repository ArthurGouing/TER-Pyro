#ifndef _FINITEVOLUME_H

#include <string>
#include "Dense"
#include "Sparse"
#include "Function.h"

class FiniteVolume
{
private:
	Eigen::SparseMatrix<double> _mat_flux;
	Eigen::VectorXd _BC_RHS;
	Function* _fct;
	DataFile* _df;
	//Adapt_maillage* _adm;
	//Vecteur solution exacte
  	Eigen::VectorXd _solexact;

public:
	// Constructeur
	FiniteVolume(Function* function, DataFile* data_file);
	void Build_flux_mat();
	void Build_BC_RHS(const double& t);
	Eigen::VectorXd ExactSolution(double t);
	const Function* Get_fct() const {return _fct;};
	const Eigen::SparseMatrix<double> Get_mat_flux() const {return _mat_flux;};
	const Eigen::VectorXd Get_BC_RHS() const {return _BC_RHS;};
};

#define _FINITEVOLUME_H
#endif
