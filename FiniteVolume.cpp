#ifndef _FINITEVOLUME_CPP

#include "FiniteVolume.h"
#include <iostream>

using namespace std;
using namespace Eigen;

// Constructeur
FiniteVolume::FiniteVolume(Function* function, DataFile* data_file ,Mesh_Adapt * adapt_mesh) :
_fct(function), _df(data_file), _adm(adapt_mesh) //il faut mettre en entré le truc p.5
{
}



//Construction de la matrice H
void FiniteVolume::Build_flux_mat(VectorXd rho, VectorXd rhostar)
{
	vector<Triplet<double>> triplets;
	int Nx=_df->Get_Nx();
	int Ny= _df->Get_Ny();
	double dt = _df->Get_dt();
	double dx = _df->Get_dx();
	VectorXd Dy =_adm->Get_Dy(), sigma(Nx*Ny);


	//Calcul de sigmak
	double lambdapv=_df->Get_lambdapv(), cpv = _df->Get_cpv();
	for (int k=0; k<Nx*Ny; ++k)
	{
		sigma(k)=lambdapv/(rho(k)*cpv);
	}


	//Calcul de Idtilde, ajout de sigma dans H
	SparseMatrix<double> Idtilde;
	Idtilde.resize(Nx*Ny,Nx*Ny);
	Idtilde.setIdentity();
	vector<Triplet<double>> tripletsId;
	for (int i=0; i<Nx*Ny; ++i)
	{
		tripletsId.push_back({i,i,rhostar(i)/rho(i)});
	}
	Idtilde.setFromTriplets(tripletsId.begin(), tripletsId.end());


	//------------------------------ Fx, Fy ----------------------------------//
	double Fx;
	VectorXd Fy(Ny) , My(Ny-1);
	Fx=1./pow(dx,2);
	for (int i=0; i<Ny-1 ; i++)
	{
		Fy(i)=1./Dy(i);
		My(i)=2./(Dy(i)+Dy(i+1));
	}
	Fy(Ny-1)=1./Dy(Ny-1);

	//------------------Id+sigma*dt*H ou Idtilde+sigmak*dt*H-------------------//
	_mat_flux.resize(Nx*Ny,Nx*Ny);

	//sigmak*H
	for (int j=1; j<=Ny; ++j)
	{
		for (int k=(j-1)*Nx+1; k<=j*Nx-1; ++k) //on se balade suivant les arrêtes suivant une ligne
		{
			triplets.push_back({k-1,k-1,sigma(k-1)*Fx});     //coef qu'il gagne de k
			triplets.push_back({k-1,k,sigma(k-1)*-Fx});     //coef qu'il perd et donne à k
			triplets.push_back({k,k,sigma(k)*Fx});        //coef qu'il gagne de k-1
			triplets.push_back({k,k-1,sigma(k)*-Fx});     //coef qu'il perd et donne à k-1
		}
	}

	for (int i=1; i<=Nx; ++i)
	{
		for (int k=1; k<=Ny-1; ++k) //on se balade suivant les arrêtes suivant une colonne
		{
			triplets.push_back({i+(k-1)*Nx-1,i+(k-1)*Nx-1,sigma(i+(k-1)*Nx-1)*Fy(k-1)*My(k-1)});
			triplets.push_back({i+(k-1)*Nx-1,i+(k-1)*Nx+Nx-1,sigma(i+(k-1)*Nx-1)*-Fy(k-1)*My(k-1)});
			triplets.push_back({i+(k-1)*Nx+Nx-1,i+(k-1)*Nx+Nx-1,sigma(i+(k-1)*Nx+Nx-1)*Fy(k)*My(k-1)});
			triplets.push_back({i+(k-1)*Nx+Nx-1,i+(k-1)*Nx-1,sigma(i+(k-1)*Nx+Nx-1)*-Fy(k)*My(k-1)});
		}
	}
	_mat_flux.setFromTriplets(triplets.begin(), triplets.end());

	//Idtilde + sigmak*dt*H
	_mat_flux=Idtilde+dt*_mat_flux;

	// cout << "-------------------------------" << endl;
	// cout << "_mat_flux (build with triplets) = " << endl;
	// cout << _mat_flux << endl;
	// cout << "-------------------------------" << endl;

}



void FiniteVolume::Build_BC_RHS(const double& t, VectorXd rho, VectorXd rhostar) //Validé
{
	int Nx=_df->Get_Nx();
	int Ny= _df->Get_Ny();
	double dt = _df->Get_dt();
	VectorXd Dy =_adm->Get_Dy();
	double lambdapv=_df->Get_lambdapv(), cpv = _df->Get_cpv();
	double Lm = _df->Get_Lm();
	double T0 = _df->Get_T0();
	_BC_RHS.resize(Nx*Ny);

	//Partie non pyro
	for (int i=0; i<=Nx*Ny-Nx-1; i++)
	{
		_BC_RHS(i)=0.;
	}
	for (int i=Nx*Ny-Nx; i<=Nx*Ny-1; i++)
	{
		_BC_RHS(i)=(dt*lambdapv/(rho(i)*cpv*Dy(Dy.size()-1)))*_fct->SourceFunction(t);
	}

	//Partie pyro
	for (int i=0; i<=Nx*Ny-1; i++)
	{
		_BC_RHS(i)+=-(1.-(rhostar(i)/rho(i)))*((Lm/cpv)+T0);
	}

	// cout << "-------------------------------" << endl;
	// cout << "_BC_RHS = " << endl;
	// cout << _BC_RHS << endl;
	// cout << "-------------------------------" << endl;
}


#define _FINITEVOLUME_CPP
#endif
