#ifndef _FINITEVOLUME_H

#include <string>
#include "Function.h"
#include "MeshAdapt.h"

class FiniteVolume
{
private:
	Eigen::SparseMatrix<double> _mat_flux;
	Eigen::VectorXd _BC_RHS;
	Function* _fct;
	DataFile* _df;
	Mesh_Adapt* _adm;
	//Vecteur solution exacte
  Eigen::VectorXd _solexact;

public:
	// Constructeur
	FiniteVolume(Function* function, DataFile* data_file, Mesh_Adapt * adapt_mesh);
	void Build_flux_mat(Eigen::VectorXd rho, Eigen::VectorXd rhostar);
	void Build_BC_RHS(const double& t, Eigen::VectorXd rho, Eigen::VectorXd rhostar);

	void Build_flux_mat_ALE(Eigen::VectorXd rho, Eigen::VectorXd rhostar);
	void Build_BC_RHS_ALE(const double& t, Eigen::VectorXd rho, Eigen::VectorXd rhostar);

	Eigen::VectorXd ExactSolution(double t);
	Function* Get_fct() const {return _fct;};
	const Eigen::SparseMatrix<double> Get_mat_flux() const {return _mat_flux;};
	const Eigen::VectorXd Get_BC_RHS() const {return _BC_RHS;};
};

#define _FINITEVOLUME_H
#endif
