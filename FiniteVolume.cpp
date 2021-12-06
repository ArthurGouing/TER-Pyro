#ifndef _FINITEVOLUME_CPP

#include "FiniteVolume.h"
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructeur
FiniteVolume::FiniteVolume(Function* function, DataFile* data_file) :
_fct(function), _df(data_file) //il faut mettre en entré le truc p.5
{
	Build_flux_mat();
}



//Construction de la matrice H
void FiniteVolume::Build_flux_mat()
{
	vector<Triplet<double>> triplets;
	int Nx=_df->Get_Nx();
	int Ny= _df->Get_Ny();
	double dt = _df->Get_dt();
	double dx = _df->Get_dx();
	VectorXd Dy=_adm->Get_Dy();
	double sigma = _df->Get_sigma(), Fx;
	VectorXd 	Fy(Ny);
	_mat_flux.resize(Nx*Ny,Nx*Ny);
	Fx=1./pow(dx,2);
	
	for (int i=0; i<Ny ; i++)
	{
		Fy(i)=1./pow(Dy(i),2);
	}

	for (int j=1; j<=Ny; ++j)
	{
		for (int k=(j-1)*Nx+1; k<=j*Nx-1; ++k) //on se balade suivant les arrêtes suivant une ligne
		{
			triplets.push_back({k-1,k-1,Fx});     //coef qu'il gagne de k
			triplets.push_back({k-1,k,-Fx});     //coef qu'il perd et donne à k
			triplets.push_back({k,k,Fx});        //coef qu'il gagne de k-1
			triplets.push_back({k,k-1,-Fx});     //coef qu'il perd et donne à k-1
		}
	}

	for (int i=1; i<=Nx; ++i)
	{
		for (int k=1; k<=Ny-1; ++k) //on se balade suivant les arrêtes suivant une colonne
		{
			triplets.push_back({i+(k-1)*Nx-1,i+(k-1)*Nx-1,Fy(k-1)});
			triplets.push_back({i+(k-1)*Nx-1,i+(k-1)*Nx+Nx-1,-Fy(k-1)});
			triplets.push_back({i+(k-1)*Nx+Nx-1,i+(k-1)*Nx+Nx-1,Fy(k)});
			triplets.push_back({i+(k-1)*Nx+Nx-1,i+(k-1)*Nx-1,-Fy(k)});
		}
	}
	_mat_flux.setFromTriplets(triplets.begin(), triplets.end());

	// cout << "-------------------------------" << endl;
	// cout << "_mat_flux (build with triplets) = " << endl;
	// cout << _mat_flux << endl;
	// cout << "-------------------------------" << endl;

}



void FiniteVolume::Build_BC_RHS(const double& t) //Validé
{
	int Nx=_df->Get_Nx();
	int Ny= _df->Get_Ny();
	double dx = _df->Get_dx();
	VectorXd Dy=_adm->Get_Dy();
	double xmin = _df->Get_xmin();
	double ymin = _df->Get_ymin();
	double ppv = _df->Get_ppv();
	double cpv = _df->Get_cpv();
	_BC_RHS.resize(Nx*Ny);

	for (int i=0; i<=Nx*Ny-Nx-1; i++)
	{
		_BC_RHS(i)=0.;
	}
	for (int i=Nx*Ny-Nx; i<=Nx*Ny-1; i++)
	{
		_BC_RHS(i)=(1./(ppv*cpv*Dy(i)))*_fct->SourceFunction(t);
	}

	// cout << "-------------------------------" << endl;
	// cout << "_BC_RHS = " << endl;
	// cout << _BC_RHS << endl;
	// cout << "-------------------------------" << endl;
}


#define _FINITEVOLUME_CPP
#endif
