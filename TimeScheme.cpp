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
}


/*void TimeScheme::ExactSolution(double t)
{
	int Nx=_df->Get_Nx();
	int Ny= _df->Get_Ny();
	double dx = _df->Get_dx();
	double dy = _df->Get_dy();
	_sol.resize(Nx*Ny);

	for (int j=0; j<Ny-1; ++j)
	{
		for (int i=0; i<Nx-1; ++i)
		{
			_solexact(i)=_fin_vol->Get_fct()->ExactSolution((i+1)*dx,(j+1)*dx,t);
		}
	}
}
*/

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


const Eigen::VectorXd & TimeScheme::GetSolution() const
{
  return _sol;
}


EulerScheme::EulerScheme(DataFile* data_file, FiniteVolume* fin_vol) :
TimeScheme(data_file,fin_vol)
{
}


void EulerScheme::Advance()
{

	double sigma=_df->Get_sigma(), dt=_df->Get_dt();
	Eigen::SparseMatrix<double> H=_fin_vol->Get_mat_flux();
	Eigen::VectorXd RHS=_fin_vol->Get_BC_RHS();

	_sol=_sol+dt*(-sigma*H*_sol+RHS);
}


ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, FiniteVolume* fin_vol) :
TimeScheme(data_file,fin_vol)
{
	std::cout << "Build time implicit scheme class." << std::endl;
  std::cout << "-------------------------------------------------" << std::endl;
  SparseMatrix<double> H=_fin_vol->Get_mat_flux();
  double dt=_df->Get_dt(), sigma=_df->Get_sigma();
  SparseMatrix<double> Id;
  Id.resize(H.rows(),H.rows());
  Id.setIdentity();

  _IdplusdtsigmaH=Id+dt*sigma*H;
  _solver_direct.analyzePattern(_IdplusdtsigmaH);
  _solver_direct.factorize(_IdplusdtsigmaH);

	// cout << "-------------------------------" << endl;
	// cout << "_IdplusdtsigmaH (build with triplets) = " << endl;
	// cout << _IdplusdtsigmaH << endl;
	// cout << "-------------------------------" << endl;
}


void ImplicitEulerScheme::Advance()
{
	double dt=_df->Get_dt();
	_t=_t+dt;
	_fin_vol->Build_BC_RHS(_t);
	Eigen::VectorXd BC_RHS=_fin_vol->Get_BC_RHS();
	Eigen::VectorXd b;
	b=_sol+dt*BC_RHS;

	_sol=_solver_direct.solve(b);

	// cout << "-------------------------------" << endl;
	// cout << "_sol = " << endl;
	// cout << _sol << endl;
	// cout << "-------------------------------" << endl;
}

#define _TIME_SCHEME_CPP
#endif
