#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>
#include <fstream>
#include <cmath>


using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, FiniteVolume* fin_vol) :
_df(data_file), _fin_vol(fin_vol), _t(_df->Get_t0()), _sol(data_file), _solold(data_file)   //whattttt ???,
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
	_sol.Get_T().resize(Nx*Ny);//inutile et ne compile pas ?????

	for (int j=0; j<Ny; ++j)
	{
		for (int i=0; i<Nx; ++i)
		{
			//_sol.T(j*Nx+i) =_fin_vol->Get_fct()->InitialCondition((i+1)*dx+xmin,(j+1)*dy+ymin);
			_sol.T(j*Nx+i) =  _fin_vol->Get_fct()->InitialCondition((i+1)*dx+xmin,(j+1)*dy+ymin);

		}
	}
	// cout << "-------------------------------" << endl;
	// cout << "_sol0 (build with triplets) = " << endl;
	// cout << _sol << endl;
	// cout << "-------------------------------" << endl;

	_sol.Get_rho().resize(Nx*Ny);//inutile

	//Ajout d'un vecteur masse volumique
	for (int j=0; j<Ny; ++j)
	{
		for (int i=0; i<Nx; ++i)
		{
			//_sol.rho(j*Nx+i)=_fin_vol->Get_fct()->InitialConditionrho((i+1)*dx+xmin,(j+1)*dy+ymin);
			_sol.rho(j*Nx+i) = _fin_vol->Get_fct()->InitialConditionrho((i+1)*dx+xmin,(j+1)*dy+ymin);
		}
	}
}


void TimeScheme::Update_Told_rhoold()
{
	_solold.T=_sol.T;
	_solold.rho=_sol.rho;
}


ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, FiniteVolume* fin_vol) :
TimeScheme(data_file,fin_vol)
{
	std::cout << "                                                  " << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "-------- Build time implicit scheme class --------" << std::endl;
  std::cout << "--------------------------------------------------" << std::endl;
}



void ImplicitEulerScheme::Advance()
{
	double dt=_df->Get_dt();
	_t=_t+dt;


	//Calcul de _rhostar
	Eigen::VectorXd Arr;
	Arr=_fin_vol->Get_fct()->Arrhenius(_sol.rho,_sol.T);
	_sol.rhostar=_sol.rho+dt*Arr; // ne compile pas ?, il faut une fonction set ???


	//Calcul de Tn+1
	_fin_vol->Build_flux_mat(_sol.rho,_sol.rhostar); //Build_flux_mat_and_BC_RHS(_t);
	_fin_vol->Build_BC_RHS(_t,_sol.rho,_sol.rhostar);
	Eigen::VectorXd BC_RHS=_fin_vol->Get_BC_RHS();
	SparseMatrix<double> A=_fin_vol->Get_mat_flux();
	Eigen::VectorXd b;
	_solver_direct.analyzePattern(A);
	_solver_direct.factorize(A);

	b=_sol.T+BC_RHS;
	_sol.T=_solver_direct.solve(b);


	//Calcul de rhon+1
	double Aref=_df->Get_Aref(), Ta=_df->Get_Ta(), rhov=_df->Get_rhov(), rhop=_df->Get_rhop();
	Arr=_fin_vol->Get_fct()->Arrhenius(_sol.rhostar,_sol.T);
	double B = rhov*Aref*dt/(rhov-rhop);
	for (int i=0; i<_sol.rho.size() ;i++)
	{
		_sol.rho(i)=(_sol.rho(i)+B*rhop*exp(-Ta/_sol.T(i)))/(1.+B*exp(-Ta/_sol.T(i)));//c'est la méthode rho(double n)
	}


	// cout << "-------------------------------" << endl;
	// cout << "_sol = " << endl;
	//cout << A << endl;
	//cout << "voilà b" << _t << endl;
	//cout << b << endl;
	// cout << "-------------------------------" << endl;

}



void ImplicitEulerScheme::Advance_ALE()
{
	double dt=_df->Get_dt();
	_t=_t+dt;


	//Calcul de _rhostar
	Eigen::VectorXd Arr;
	Arr=_fin_vol->Get_fct()->Arrhenius(_solold.rho,_solold.T);
	_solold.rhostar=_solold.rho+dt*Arr; // ne compile pas ?, il faut une fonction set ???

	//cout << "après _rhostar" << endl;
	//Calcul de Tn+1
	_fin_vol->Build_flux_mat_ALE(_solold); //Build_flux_mat_and_BC_RHS(_t);
	_fin_vol->Build_BC_RHS_ALE(_t,_solold);   
	Eigen::VectorXd BC_RHS=_fin_vol->Get_BC_RHS();
	SparseMatrix<double> A=_fin_vol->Get_mat_flux();
	Eigen::VectorXd b;
	_solver_direct.analyzePattern(A);
	_solver_direct.factorize(A);

	b=_solold.T+BC_RHS;
	_sol.T=_solver_direct.solve(b);

	//cout << "après Tn+1" << endl;
	//Calcul de rhon+1
	double Aref=_df->Get_Aref(), Ta=_df->Get_Ta(), rhov=_df->Get_rhov(), rhop=_df->Get_rhop();
	Arr=_fin_vol->Get_fct()->Arrhenius(_solold.rhostar,_sol.T);
	double B = rhov*Aref*dt/(rhov-rhop);
	for (int i=0; i<_sol.rho.size() ;i++)
	{
		_sol.rho(i)=(_solold.rho(i)+B*rhop*exp(-Ta/_sol.T(i)))/(1.+B*exp(-Ta/_sol.T(i)));//c'est la méthode rho(double n)
	}

	//cout << "après rhon+1" << endl;
	// cout << "-------------------------------" << endl;
	// cout << "_sol = " << endl;
	//cout << A << endl;
	//cout << "voilà b" << _t << endl;
	//cout << b << endl;
	// cout << "-------------------------------" << endl;

}



void TimeScheme::SaveSol(Solution sol, string n_sol, int n)
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
			solution << sol.T(i+j*Nx) << " ";
		}
		solution << endl;
	}
	solution.close();
}

void TimeScheme::Save_rho(Eigen::VectorXd rho , double t , std::string name_file)
{
	string n_file = name_file + to_string(t) + ".txt";
	int Nx(_df->Get_Nx()), Ny(_df->Get_Ny());
	double xmin(_df->Get_xmin()), ymin(_df->Get_ymin());
	double dx(_df->Get_dx()), dy(_df->Get_dy());
	ofstream solution_rho;
	solution_rho.open(name_file, ios::out);
	for (int i=0 ; i<Ny ;i++)
	{
		solution_rho << i*dy << rho(i*Nx) <<endl;
	}
	solution_rho.close();
}

#define _TIME_SCHEME_CPP
#endif
