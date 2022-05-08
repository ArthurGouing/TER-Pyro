#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>
#include <fstream>
#include <cmath>


using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, FiniteVolume* fin_vol, Mesh_Adapt* adm) :
_df(data_file), _fin_vol(fin_vol), _t(_df->Get_t0()), _sol(data_file), _solold(data_file), _adm(adm)  //whattttt ???,
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
	// _sol.Get_T().resize(Nx*Ny);//inutile et ne compile pas ?????
	// _solold.Get_T().resize(Nx*Ny);
	_sol.T.resize(Nx*Ny);//inutile et ne compile pas ?????
	_solold.T.resize(Nx*Ny);

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

	_sol.rho.resize(Nx*Ny);//inutile
	_solold.rho.resize(Nx*Ny);
	_sol.rhostar.resize(Nx*Ny);//inutile
	_solold.rhostar.resize(Nx*Ny);

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
	// for ( int i=0; i<_solold.rhostar.size() ;i++)
	// {
	// 	//cout << "i=" << i << "Update rhoold(i)=" << _solold.rho(i) << endl;
	// }
	// for ( int i=0; i<_solold.rhostar.size() ;i++)
	// {
	// 	cout << "i=" << i << "Update rho(i)=" << _sol.rho(i) << endl;
	// }
	_solold.rho=_sol.rho;
}


ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, FiniteVolume* fin_vol, Mesh_Adapt* adm) :
TimeScheme(data_file,fin_vol, adm)
{
	std::cout << "                                                  " << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "-------- Build time implicit scheme class --------" << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;
}



void ImplicitEulerScheme::Advance(double tn)
{
	double dt=_df->Get_dt();
	_t=tn;


	//Calcul de _rhostar
	Eigen::VectorXd Arr;
	Arr=_fin_vol->Get_fct()->Arrhenius(_sol.rho,_sol.T);
	cout << "Debut adance normale" << endl;
	for ( int i=0; i<_solold.rhostar.size() ;i++)
	{
		//cout << "i=" << i << "Advance normale rho(i)=" << _sol.rho(i) << endl;
	}
	_sol.rhostar=_sol.rho+dt*Arr; // ne compile pas ?, il faut une fonction set ???
	cout << " adance normale rhostar" << endl;
	for ( int i=0; i<_solold.rhostar.size() ;i++)
	{
		//cout << "i=" << i << "Advance normale rhostar(i)=" << _sol.rhostar(i) << endl;
	}

	//Calcul de Tn+1
	_fin_vol->Build_flux_mat(_sol.rho,_sol.rhostar); //Build_flux_mat_and_BC_RHS(_t);
	_fin_vol->Build_BC_RHS(_t,_sol.rho,_sol.rhostar);
	Eigen::VectorXd BC_RHS=_fin_vol->Get_BC_RHS();
	SparseMatrix<double> A=_fin_vol->Get_mat_flux();
	Eigen::VectorXd b;
	_solver_direct.analyzePattern(A);
	_solver_direct.factorize(A);

	for(int i=0; i<BC_RHS.size();i++)
	{
	// 	for(int j=0; j<BC_RHS.size();j++){
	// 		cout <<"A" << i<<", "<<j<<": " << A.coeff(i,j) << " " << A_n.coeff(i,j)<<endl;
	// 		// moy_1+=abs(A.coeff(i,j)-A_n.coeff(i,j));
	// 	}
	 		cout << "normalB" <<i<<" : "<< BC_RHS(i) << " " <<endl;
	// 		// moy_2+=abs(BC_RHS(i)-BC_RHS_n(i));
	 	}
	//////////////////////////////DEBUG
	// for(int i=0; i<BC_RHS.size();i++)
	// {
	// 	for(int j=0; j<BC_RHS.size();j++){
	// 		cout << "A" << i<<", "<<j<<": " << A.coeff(i,j) <<endl;
	// 		// moy_1+=abs(A.coeff(i,j)-A_n.coeff(i,j));
	// 	}
	// 		cout << "B" <<i<<" : "<< BC_RHS(i) <<endl;
	// 		// moy_2+=abs(BC_RHS(i)-BC_RHS_n(i));
	// 	}


	b=_sol.T+BC_RHS;
	_sol.T=_solver_direct.solve(b);
	for (int i=0; i<_sol.rho.size() ;i++)
	{
		cout << "i=" << i << " Tn+1(i)=" << _sol.T(i) << endl;
	}

	//Calcul de rhon+1
	double Aref=_df->Get_Aref(), Ta=_df->Get_Ta(), rhov=_df->Get_rhov(), rhop=_df->Get_rhop();
	Arr=_fin_vol->Get_fct()->Arrhenius(_sol.rhostar,_sol.T);
	double B = rhov*Aref*dt/(rhov-rhop);
	for (int i=0; i<_sol.rho.size() ;i++)
	{
		_sol.rho(i)=(_sol.rho(i)+B*rhop*exp(-Ta/_sol.T(i)))/(1.+B*exp(-Ta/_sol.T(i)));
	}
	//cout << " adance normale rho n+1" << endl;
	for ( int i=0; i<_solold.rhostar.size() ;i++)
	{
		cout << "i=" << i << "Advance normale rho(i)=" << _sol.rho(i) << endl;
	}


	// cout << "-------------------------------" << endl;
	// cout << "_sol = " << endl;
	//cout << A << endl;
	//cout << "voilà b" << _t << endl;
	//cout << b << endl;
	// cout << "-------------------------------" << endl;

}



void ImplicitEulerScheme::Advance_ALE(double tn)
{
	double dt=_df->Get_dt();
	_t=tn;
	//Calcul de _rhostar
	int Nx = _df->Get_Nx();
	Eigen::VectorXd Arr  = _fin_vol->Get_fct()->Arrhenius(_solold.rho,_solold.T);;
	Eigen::VectorXd Dy   = _adm->Get_Dy();
	Eigen::VectorXd Dyold= _adm->Get_Dyold();
	//Debugbegin
	for ( int i=0; i<Dyold.size() ;i++)
	{
		//cout << "Dyold("<<i<<")="<<Dyold(i) << endl;
	}
	for ( int i=0; i<Dyold.size() ;i++)
	{
		//cout << "Dy("<<i<<")="<<Dy(i) << endl;
	}
	//Debugend
	for ( int i=0; i<_solold.rhostar.size() ;i++)
	{
		//cout << "i=" << i << " rho(i)=" << _solold.rho(i) << endl;
		_solold.rhostar(i)=_solold.rho(i)*(Dyold(i/Nx)/Dy(i/Nx)) + (dt/Dy(i/Nx))*(advecrho()(i)+Arr(i)*Dyold(i/Nx)); // On chance F par le terme d'advection
		//cout << "i=" << i << " rhostar(i)=" << _solold.rhostar(i) << endl;
	}

	cout << "après _rhostarssss" << endl;
	//Calcul de Tn+1
	_fin_vol->Build_flux_mat_ALE(_solold); //Build_flux_mat_and_BC_RHS(_t);
	cout << "apres Build_flux_mat_ALE" << endl;
	_fin_vol->Build_BC_RHS_ALE(_t,_solold);
	cout << "avant get" << endl;
	Eigen::VectorXd BC_RHS=_fin_vol->Get_BC_RHS();
	SparseMatrix<double> A=_fin_vol->Get_mat_flux();
	cout << "apres get" << endl;

/***************************** Pour deguguer **********************************/
	// _fin_vol->Build_flux_mat(_solold.rho,_solold.rhostar);
	// _fin_vol->Build_BC_RHS(_t,_solold.rho,_solold.rhostar);
	// Eigen::VectorXd BC_RHS_n=_fin_vol->Get_BC_RHS();
	// SparseMatrix<double> A_n=_fin_vol->Get_mat_flux();
	// float moy_1=0;
	// float moy_2=0;
	for(int i=0; i<BC_RHS.size();i++)
	{
	// 	for(int j=0; j<BC_RHS.size();j++){
	// 		cout <<"A" << i<<", "<<j<<": " << A.coeff(i,j) << " " << A_n.coeff(i,j)<<endl;
	// 		// moy_1+=abs(A.coeff(i,j)-A_n.coeff(i,j));
	// 	}
	 		cout << "ALEB" <<i<<" : "<< BC_RHS(i) << " " <<endl;
	// 		// moy_2+=abs(BC_RHS(i)-BC_RHS_n(i));
	 	}
		// cout << "Erreur mmoyenne de A : "<<moy_1/(BC_RHS.size()*BC_RHS.size())<<endl;
		// cout << "Erreur mmoyenne de BC_RHS : "<<moy_2/(BC_RHS.size())<<endl;
/*
Il y a une autre erreur dans le TimeScheme car on a le même résultat si on prend
les matrices ALE ou les matrices normales dans le Advance_ALE, peut être juste
une mauvaise initialisation de sol old.
Donc les matrices ALE sont correct (peut etre pas les temres UpWind)
************************* Fin débuguage ***************************************/

		Eigen::VectorXd b;
		_solver_direct.analyzePattern(A);
		_solver_direct.factorize(A);
		cout << "avant solver" << endl;
		b=_solold.T+BC_RHS;
		_sol.T=_solver_direct.solve(b);
		for (int i=0; i<_sol.rho.size() ;i++)
		{
			cout << "i=" << i << " Tn+1advec(i)=" << _sol.T(i) << endl;
		}
		cout << "après Tn+1" << endl;
		//Calcul de rhon+1

		double Aref=_df->Get_Aref(), Ta=_df->Get_Ta(), rhov=_df->Get_rhov(), rhop=_df->Get_rhop();
		Arr=_fin_vol->Get_fct()->Arrhenius(_solold.rhostar,_sol.T);
		double B = rhov*Aref*dt/(rhov-rhop);
		for (int i=0; i<_sol.rho.size() ;i++)
		{
			_sol.rho(i)=((Dyold(i/Nx)/Dy(i/Nx))*_solold.rho(i)+(dt/Dy(i/Nx))*advecrho()(i)+B*rhop*exp(-Ta/_sol.T(i)))/(1.+B*exp(-Ta/_sol.T(i)));//c'est la méthode rho(double n)
			cout << "i=" << i << " ALE rhon+1(i)=" << _sol.rho(i) << endl;
		}

		//cout << "après rhon+1" << endl;
		// cout << "-------------------------------" << endl;
		// cout << "_sol = " << endl;
		//cout << A << endl;
		//cout << "voilà b" << _t << endl;
		//cout << b << endl;
		// cout << "-------------------------------" << endl;

	}


	VectorXd TimeScheme::advecrho()
	{
		//aide
		int Nx=_df->Get_Nx();
		int Ny= _df->Get_Ny();
		VectorXd advec(Nx*Ny), c; //vitesse d'advection
		c=_adm->Get_vitesse();

		for (int i=0; i<=Nx-1; i++) //Première ligne
		{
			//Terme d'UpWind
			//cout << "c("<<(i/Nx)+1<<")="<<c((i/Nx)+1) << endl;
			advec(i)= max(0,c((i/Nx)+1)) * _solold.rho(i)
			- max(0,-c((i/Nx)+1)) * _solold.rho(i+Nx);
			//cout << "dans advec premiere ligne " << "advec("<<i<<")="<< advec(i) << endl;
		}
		for (int i=Nx; i<=Nx*Ny-Nx-1; i++) //Lignes intermédiaires
		{
			//Terme d'UpWind
			//cout << "c("<<(i/Nx)+1<<")="<<c((i/Nx)+1) << endl;
			//cout << "c("<<(i/Nx)<<")="<<c((i/Nx)) << endl;
			advec(i)=max(0,c((i/Nx)+1)) * _solold.rho(i)
			- max(0,-c((i/Nx)+1)) * _solold.rho(i+Nx)
			+ max(0,-c(i/Nx)) * _solold.rho(i)
			- max(0,c(i/Nx)) * _solold.rho(i-Nx);
		}
		for (int i=Nx*Ny-Nx; i<=Nx*Ny-1; i++) //Dernière ligne
		{
			//Terme d'UpWind
			//cout << "c("<<(i/Nx)<<")="<<c((i/Nx)) << endl;
			advec(i)= max(0,-c(i/Nx)) * _solold.rho(i)
			- max(0,c(i/Nx)) * _solold.rho(i-Nx);
		}
		return advec;
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

	double TimeScheme::max(double a, double b)
	{
		if (a>b)
		{
			return a;
		}
		else
		{
			return b;
		}
	}

	#define _TIME_SCHEME_CPP
	#endif
