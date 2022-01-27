#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>
#include <fstream>
#include <cmath>


using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, FiniteVolume* fin_vol) :
_df(data_file), _fin_vol(fin_vol), _t(_df->Get_t0())
{
	InitialCondition();
}


// Destructeur
TimeScheme::~TimeScheme()
{
}


void TimeScheme::InitialCondition()
{
	int Nx=_df->Get_Nx();
	int Ny= _df->Get_Ny();
	double dx = _df->Get_dx();
	double dy = _df->Get_dy();
	double xmin = _df->Get_xmin();
  double ymin = _df->Get_ymin();
	_sol.resize(Nx*Ny);

	for (int j=0; j<Ny; ++j)
	{
		for (int i=0; i<Nx; ++i)
		{
			_sol(j*Nx+i)=_fin_vol->Get_fct()->InitialCondition((i+1)*dx+xmin,(j+1)*dy+ymin);
		}
	}
	// cout << "-------------------------------" << endl;
	// cout << "_sol0 (build with triplets) = " << endl;
	// cout << _sol << endl;
	// cout << "-------------------------------" << endl;

	_rho.resize(Nx*Ny);

	//Ajout d'un vecteur masse volumique
	for (int j=0; j<Ny; ++j)
	{
		for (int i=0; i<Nx; ++i)
		{
			_rho(j*Nx+i)=_fin_vol->Get_fct()->InitialConditionrho((i+1)*dx+xmin,(j+1)*dy+ymin);
		}
	}
}


ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, FiniteVolume* fin_vol) :
TimeScheme(data_file,fin_vol)
{
	std::cout << "Build time implicit scheme class." << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
}


VectorXd ImplicitEulerScheme::rhostarexp(VectorXd rho, VectorXd sol)
{
	Eigen::VectorXd Arr, rhonpun;
	double dt = _df->Get_dt();

	Arr=_fin_vol->Get_fct()->Arrhenius(rho,sol);
	rhonpun=rho+dt*Arr;

	return rhonpun;
}


void ImplicitEulerScheme::Advance()
{
	double dt=_df->Get_dt();
	_t=_t+dt;


	//Calcul de _rhostar
	_rhostar=rhostarexp(_rho, _sol);


	//Calcul de Tn+1
	_fin_vol->Build_flux_mat(_rho,_rhostar); //Build_flux_mat_and_BC_RHS(_t);
	_fin_vol->Build_BC_RHS(_t,_rho,_rhostar);
	Eigen::VectorXd BC_RHS=_fin_vol->Get_BC_RHS();
	SparseMatrix<double> A=_fin_vol->Get_mat_flux();
	Eigen::VectorXd b;
	_solver_direct.analyzePattern(A);
	_solver_direct.factorize(A);

	b=_sol+BC_RHS;
	_sol=_solver_direct.solve(b);


	//Calcul de rhon+1
	_rho=rhostarexp(_rhostar, _sol);

	// cout << "-------------------------------" << endl;
	// cout << "_sol = " << endl;
	//cout << A << endl;
	//cout << "voilà b" << _t << endl;
	//cout << b << endl;
	// cout << "-------------------------------" << endl;
}


const Eigen::VectorXd & TimeScheme::GetSolution() const
{
  return _sol;
}


void TimeScheme::SaveSol(Eigen::VectorXd sol, string n_sol, int n)
{
	string n_file = _df->Get_results() + "/" + n_sol + to_string(n) + ".vtk";
	ofstream solution;
	solution.open(n_file, ios::out);
	int Nx(_df->Get_Nx()), Ny(_df->Get_Ny());
	double xmin(_df->Get_xmin()), ymin(_df->Get_ymin());
	double dx(_df->Get_dx()), dy(_df->Get_dy());

	solution << "# vtk DataFile Version 3.0" << endl;
	solution << "sol" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET STRUCTURED_POINTS" << endl;
	solution << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << endl;
	solution << "ORIGIN " << xmin << " " << ymin << " " << 0 << endl;
	solution << "SPACING " << dx << " " << dy << " " << 1 << endl;;
	solution << "POINT_DATA " << Nx*Ny << endl;
	solution << "SCALARS sol float" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	for(int j=0; j<Ny; ++j)
	{
		for(int i=0; i<Nx; ++i)
		{
			solution << sol(i+j*Nx) << " ";
		}
		solution << endl;
	}
	solution.close();
}

#define _TIME_SCHEME_CPP
#endif
