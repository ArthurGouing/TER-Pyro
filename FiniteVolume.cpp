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
	VectorXd Dy =_adm->Get_Dyold(), sigma(Nx*Ny); ////ON a mis Dyold



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

	cout << "-------------------------------" << endl;
	cout << "_mat_flux advec de base (build with triplets) = " << endl;
	cout << _mat_flux << endl;
	cout << "-------------------------------" << endl;

}




void FiniteVolume::Build_BC_RHS(const double& t, VectorXd rho, VectorXd rhostar) //Validé
{
	int Nx=_df->Get_Nx();
	int Ny= _df->Get_Ny();
	double dt = _df->Get_dt();
	VectorXd Dy =_adm->Get_Dyold(); ///// On a mis Dyold !!!!!!!!!!!!!
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
		//_BC_RHS(i)=0.0;
	//	cout << "Debug second membre partie non pyro " << "B" << i << " "<< dt << " " << lambdapv << " "  << rho(i) << " "  << cpv<< " "  <<Dy(Dy.size()-1) << " "  << _fct->SourceFunction(t) << endl;
	}

	//Partie pyro
	for (int i=0; i<=Nx*Ny-1; i++)
	{
		//_BC_RHS(i)+=(1.-(rhostar(i)/rho(i)))*((Lm/cpv)-T0); ///////////!!!!!!
		_BC_RHS(i)+=-(1.-(rhostar(i)/rho(i)))*((Lm/cpv)+T0);
		//cout << "Debug second membre partie pyro " << "B" << i << " "<< _BC_RHS(i) << endl;
	}

	// cout << "-------------------------------" << endl;
	// cout << "_BC_RHS = " << endl;
	// cout << _BC_RHS << endl;
	// cout << "-------------------------------" << endl;
}


void FiniteVolume::Build_flux_mat_ALE(Solution sol)
{
	vector<Triplet<double>> triplets;
	int Nx    =_df->Get_Nx();
	int Ny    = _df->Get_Ny();
	double dt = _df->Get_dt();
	double dx = _df->Get_dx();
	VectorXd Dy     =_adm->Get_Dy(), sigma(Nx*Ny);
	VectorXd Dystar =_adm->Get_Dystar();
	VectorXd Dyold  =_adm->Get_Dyold();


	//Calcul de sigmak
	cout << "avant sigmak" << endl;
	double lambdapv=_df->Get_lambdapv(), cpv = _df->Get_cpv();
	for (int k=0; k<Nx*Ny; ++k)
	{
		sigma(k)=lambdapv/(sol.rho(k)*cpv);
	}
	cout << "après sigmak" << endl;


	//Calcul de Idtilde, ajout de sigma dans H
	SparseMatrix<double> Idtilde;
	Idtilde.resize(Nx*Ny,Nx*Ny);
	Idtilde.setIdentity();
	vector<Triplet<double>> tripletsId;
	for (int i=0; i<Nx*Ny; ++i)
	{
		//cout << "i=" <<i << endl;
		//cout <<  sol.rho(i)  << " " << sol.rhostar(i)<< endl;
		tripletsId.push_back({i,i,(sol.rhostar(i)/sol.rho(i))*(Dy(i/Nx)/Dyold(i/Nx))}); // Verif le i/Dx
	}
	Idtilde.setFromTriplets(tripletsId.begin(), tripletsId.end());

	cout << "apres Idtilde" << endl;
	//------------------------------ Fx, Fy ----------------------------------//
	double Fx;
	VectorXd Fy(Ny);
	VectorXd Fystar(Ny);
	Fx=1./pow(dx,2);
	for (int i=0; i<Ny ; i++)
	{
		Fy(i)=1./pow(Dy(i),2);
		Fystar(i)=1./pow(Dystar(i),2);
	}
	//------------------Id+sigma*dt*H ou Idtilde+sigmak*dt*H-------------------//
	_mat_flux.resize(Nx*Ny,Nx*Ny);

	//sigmak*H
	cout << "avant sigmak*H" << endl;
	for (int j=1; j<=Ny; ++j)
	{
		for (int k=(j-1)*Nx+1; k<=j*Nx-1; ++k) //on se balade suivant les arrêtes suivant une ligne
		{

			triplets.push_back({k-1,k-1,sigma(k-1)*Fx * Dystar((k-1)/Nx)/Dyold((k-1)/Nx)});     //coef qu'il gagne de k
			triplets.push_back({k-1,k,  sigma(k-1)*-Fx* Dystar((k-1)/Nx)/Dyold((k-1)/Nx)});     //coef qu'il perd et donne à k
			triplets.push_back({k,  k,  sigma(k)*Fx   * Dystar(k/Nx)/Dyold(k/Nx)});              //coef qu'il gagne de k-1
			triplets.push_back({k,  k-1,sigma(k)*-Fx  * Dystar(k/Nx)/Dyold(k/Nx)});              //coef qu'il perd et donne à k-1
		}
	}
	cout << "avant sigmak*H2" << endl;
	for (int i=1; i<=Nx; ++i)
	{
		for (int k=1; k<=Ny-1; ++k) //on se balade suivant les arrêtes suivant une colonne
		{
			// triplets.push_back({i+(k-1)*Nx-1,i+(k-1)*Nx-1 ,sigma(i+(k-1)*Nx-1)* 1./(Dyold(k-1)*Dystar(k-1)) });
			// triplets.push_back({i+(k-1)*Nx-1,i+(k)*Nx-1   ,sigma(i+(k-1)*Nx-1)*-1./(Dyold(k-1)*Dystar(k-1)) });
			// triplets.push_back({i+(k)*Nx-1,  i+(k)*Nx-1   ,sigma(i+(k-1)*Nx+Nx-1)* 1./(Dyold(k)*Dystar(k-1)) }); //k-1 au lieu de k ??
			// triplets.push_back({i+(k)*Nx-1,  i+(k-1)*Nx-1 ,sigma(i+(k-1)*Nx+Nx-1)*-1./(Dyold(k)*Dystar(k-1)) }); //k-1 au lieu de k ?? et cest pas dystar

			// triplets.push_back({i+(k-1)*Nx-1,i+(k-1)*Nx-1 ,sigma(i+(k-1)*Nx-1)* 1./(Dyold(k-1)*((Dystar(k-1)+Dystar(k))/2)) });
			// triplets.push_back({i+(k-1)*Nx-1,i+(k)*Nx-1   ,sigma(i+(k-1)*Nx-1)*-1./(Dyold(k-1)*((Dystar(k-1)+Dystar(k))/2)) });
			// triplets.push_back({i+(k)*Nx-1,  i+(k)*Nx-1   ,sigma(i+(k-1)*Nx+Nx-1)* 1./(Dyold(k)*((Dystar(k-1)+Dystar(k))/2)) }); //k-1 au lieu de k ??
			// triplets.push_back({i+(k)*Nx-1,  i+(k-1)*Nx-1 ,sigma(i+(k-1)*Nx+Nx-1)*-1./(Dyold(k)*((Dystar(k-1)+Dystar(k))/2)) });
			triplets.push_back({i+(k-1)*Nx-1,i+(k-1)*Nx-1 ,sigma(i+(k-1)*Nx-1)* 1./(Dyold((i+(k-1)*Nx-1)/Nx)*((Dystar(k-1)+Dystar(k))/2)) });
			triplets.push_back({i+(k-1)*Nx-1,i+(k)*Nx-1   ,sigma(i+(k-1)*Nx-1)*-1./(Dyold((i+(k-1)*Nx-1)/Nx)*((Dystar(k-1)+Dystar(k))/2)) });
			triplets.push_back({i+(k)*Nx-1,  i+(k)*Nx-1   ,sigma(i+(k-1)*Nx+Nx-1)* 1./(Dyold((i+(k)*Nx-1)/Nx)*((Dystar(k-1)+Dystar(k))/2)) }); //k-1 au lieu de k ??
			triplets.push_back({i+(k)*Nx-1,  i+(k-1)*Nx-1 ,sigma(i+(k-1)*Nx+Nx-1)*-1./(Dyold((i+(k)*Nx-1)/Nx)*((Dystar(k-1)+Dystar(k))/2)) });

		}
	}
	cout << "avant sigmak*H3" << endl;
	_mat_flux.setFromTriplets(triplets.begin(), triplets.end());
	double som;
	for(int i=0; i<Nx*Ny;i++)
	{
		som=0.0;
		for(int j=0; j<Nx*Ny;j++)
		{
			som+=_mat_flux.coeff(i,j);
		}
		cout << "somme sur " << i <<  "=" << som << endl;
	}

	// cout << "-------------------------------" << endl;
	// cout << "_mat_flux ALE(build with triplets) = " << endl;
	// cout << _mat_flux << endl;
	// cout << "-------------------------------" << endl;
	//Idtilde + sigmak*dt*H
	_mat_flux=Idtilde+dt*_mat_flux;

	cout << "-------------------------------" << endl;
	cout << "_mat_flux ALE(build with triplets) = " << endl;
	cout << _mat_flux << endl;
	cout << "-------------------------------" << endl;

}

// ALE pour arbitraty Lagrangian-Eulerian Methods
void FiniteVolume::Build_BC_RHS_ALE(const double& t, Solution sol)
{
	int Nx=_df->Get_Nx();
	int Ny= _df->Get_Ny();
	double dt = _df->Get_dt();
	VectorXd Dy =_adm->Get_Dy();
	VectorXd Dyold  =_adm->Get_Dyold();
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
		_BC_RHS(i)=(dt*lambdapv/(sol.rho(i)*cpv*Dyold(Dy.size()-1)))*_fct->SourceFunction(t);
		//_BC_RHS(i)=0.0;
	}

	 ///a enlever
	//Partie pyro
	for (int i=0; i<=Nx*Ny-1; i++)
	{
		// _BC_RHS(i)+= (1.-(sol.rhostar(i)/sol.rho(i)) * Dy(i/Nx)/Dyold(i/Nx) )*((Lm/cpv)-T0); //////!!!!!!!!!!!!
		 _BC_RHS(i)+= -(1.-(sol.rhostar(i)/sol.rho(i)) * (Dy(i/Nx)/Dyold(i/Nx)) )*((Lm/cpv)+T0); //////!!!!!!!!!!!!
		 // if (abs(_BC_RHS(i))<=0.000000000001)
		 // {
			//  _BC_RHS(i)=0.0;
		 // }
			cout << i << "------------------------------------ " << "Dy" <<(1.-(sol.rhostar(i)/sol.rho(i)) * (Dy(i/Nx)/Dyold(i/Nx)) )*1000000
			 << "rhostar(i)"<< sol.rhostar(i) << "rho(i)" <<sol.rho(i) << "BC(i)" << _BC_RHS(i) << endl;
		//_BC_RHS(i)+= 0.0; //////!!!!!!!!!!!!
		//cout << "Debug second membre partie pyro " << "B" << i << " "<< _BC_RHS(i) << endl;
	}
	//_BC_RHS(1)=0.0;


	//Termes upwind
	VectorXd c; //vitesse d'advection
	c=_adm->Get_vitesse();
	for (int i=0; i<=Nx-1; i++) //Première ligne
	{
		_BC_RHS(i)+= max(0.0,c((i/Nx)+1)) * ((sol.T(i)-T0)-Lm/cpv)  * dt/Dyold(i/Nx)
		- max(0.0,-c((i/Nx)+1)) * ((sol.T(i+Nx)-T0)-Lm/cpv) * dt/Dyold(i/Nx) *(sol.rho(i+Nx)/sol.rho(i));

	}
	for (int i=Nx; i<=Nx*Ny-Nx-1; i++) //Lignes intermédiaires
	{
		_BC_RHS(i)+=max(0.0,c((i/Nx)+1)) * ((sol.T(i)-T0)-Lm/cpv)  * dt/Dyold(i/Nx)
		- max(0.0,-c((i/Nx)+1)) * ((sol.T(i+Nx)-T0)-Lm/cpv) * dt/Dyold(i/Nx) *(sol.rho(i+Nx)/sol.rho(i))
		+ max(0.0,-c(i/Nx)) * ((sol.T(i)-T0)-Lm/cpv)   * dt/Dyold(i/Nx)
	  - max(0.0,c(i/Nx)) * ((sol.T(i-Nx)-T0)-Lm/cpv) * dt/Dyold(i/Nx) *(sol.rho(i-Nx)/sol.rho(i));
cout <<		max(0.0,c((i/Nx)+1)) << " " << ((sol.T(i)-T0)-Lm/cpv) <<  " " << dt/Dyold(i/Nx)  << endl;;
cout <<		- max(0.0,-c((i/Nx)+1)) << " " << ((sol.T(i+Nx)-T0)-Lm/cpv) <<  " " << dt/Dyold(i/Nx)  << " " << ( sol.rho(i+Nx)/sol.rho(i)) << endl;;
cout <<		max(0.0,-c((i/Nx))) << " " << ((sol.T(i)-T0)-Lm/cpv) <<  " " << dt/Dyold(i/Nx)  << endl;;
cout <<		- max(0.0,c((i/Nx))) << " " << ((sol.T(i-Nx)-T0)-Lm/cpv) <<  " " << dt/Dyold(i/Nx)  << " " << ( sol.rho(i-Nx)/sol.rho(i)) << endl;;

 // if (abs(_BC_RHS(i))<=0.0000001)
	// {
	// 	_BC_RHS(i)=0.0;
	// }
	cout << "ci " << max(0,c((i/Nx)+1)) << "ci+1 " << max(0,-c((i/Nx)+1)) << endl;
}
	for (int i=Nx*Ny-Nx; i<=Nx*Ny-1; i++) //Dernière ligne
	{
		_BC_RHS(i)+= max(0.0,-c(i/Nx)) * ((sol.T(i)-T0)-Lm/cpv)   * dt/Dyold(i/Nx)
		- max(0.0,c(i/Nx)) * ((sol.T(i-Nx)-T0)-Lm/cpv) * dt/Dyold(i/Nx) *(sol.rho(i-Nx)/sol.rho(i));
	}
//
// VectorXd c; //vitesse d'advection
// c=_adm->Get_vitesse();
// for (int i=0; i<=Nx-1; i++) //Première ligne
// {
// 	if (abs(Dy(i)-Dyold(i))<=0.00000000001)
// 	{
// 		_BC_RHS(i)+=0.0;
// 	}
// 	else
// 	{
// 			if (c((i/Nx)+1)>=0.00000000001)
// 	{
// 		_BC_RHS(i)+= max(0.0,c((i/Nx)+1)) * ((sol.T(i)-T0)-Lm/cpv)  * dt/Dyold(i/Nx);
// 	}
// 	if (-c((i/Nx)+1)>=0.00000000001)
// 	{
// 		_BC_RHS(i)+= 	- max(0.0,-c((i/Nx)+1)) * ((sol.T(i+Nx)-T0)-Lm/cpv) * dt/Dyold(i/Nx) *(sol.rho(i+Nx)/sol.rho(i));
// 	}
// }
// }
// for (int i=Nx; i<=Nx*Ny-Nx-1; i++) //Lignes intermédiaires
// {
// 	if (abs(Dy(i)-Dyold(i))<=0.00000000001)
// 	{
// 		_BC_RHS(i)+=0.0;
// 	}
// 	else
// 	{
// 	if (c((i/Nx)+1)>=0.000001)
// 	{
// 		_BC_RHS(i)+= max(0.0,c((i/Nx)+1)) * ((sol.T(i)-T0)-Lm/cpv)  * dt/Dyold(i/Nx);
// 	}
// 	if (-c((i/Nx)+1)>=0.000001)
// 	{
// 		_BC_RHS(i)+= 	- max(0.0,-c((i/Nx)+1)) * ((sol.T(i+Nx)-T0)-Lm/cpv) * dt/Dyold(i/Nx) *(sol.rho(i+Nx)/sol.rho(i));
// 	}
// 	if (-c((i/Nx))>=0.000001)
// 	{
// 		_BC_RHS(i)+= max(0.0,-c(i/Nx)) * ((sol.T(i)-T0)-Lm/cpv)   * dt/Dyold(i/Nx);
// 	}
// 	if (c((i/Nx))>=0.000001)
// 	{
// 		_BC_RHS(i)+= - max(0.0,c(i/Nx)) * ((sol.T(i-Nx)-T0)-Lm/cpv) * dt/Dyold(i/Nx) *(sol.rho(i-Nx)/sol.rho(i));
// 	}
// }
// cout <<		max(0.0,c((i/Nx)+1)) << " " << ((sol.T(i)-T0)-Lm/cpv) <<  " " << dt/Dyold(i/Nx)  << endl;;
// cout <<		- max(0.0,-c((i/Nx)+1)) << " " << ((sol.T(i+Nx)-T0)-Lm/cpv) <<  " " << dt/Dyold(i/Nx)  << " " << ( sol.rho(i+Nx)/sol.rho(i)) << endl;;
// cout <<		max(0.0,-c((i/Nx))) << " " << ((sol.T(i)-T0)-Lm/cpv) <<  " " << dt/Dyold(i/Nx)  << endl;;
// cout <<		- max(0.0,c((i/Nx))) << " " << ((sol.T(i-Nx)-T0)-Lm/cpv) <<  " " << dt/Dyold(i/Nx)  << " " << ( sol.rho(i-Nx)/sol.rho(i)) << endl;;
//
// // if (abs(_BC_RHS(i))<=0.0000001)
// // {
// // 	_BC_RHS(i)=0.0;
// // }
// cout << "ci " << max(0,c((i/Nx)+1)) << "ci+1 " << max(0,-c((i/Nx)+1)) << endl;
// }
// for (int i=Nx*Ny-Nx; i<=Nx*Ny-1; i++) //Dernière ligne
// {
// 	if (abs(Dy(i)-Dyold(i))<=0.00000000001)
// 	{
// 		_BC_RHS(i)+=0.0;
// 	}
// 	else
// 	{
// 	if (-c((i/Nx))>=0.00000000001)
// 	{
// 		_BC_RHS(i)+= max(0.0,-c(i/Nx)) * ((sol.T(i)-T0)-Lm/cpv)   * dt/Dyold(i/Nx);
// 	}
// 	if (c((i/Nx))>=0.00000000001)
// 	{
// 		_BC_RHS(i)+= - max(0.0,c(i/Nx)) * ((sol.T(i-Nx)-T0)-Lm/cpv) * dt/Dyold(i/Nx) *(sol.rho(i-Nx)/sol.rho(i));
// 	}
// }
//}

	// cout << "-------------------------------" << endl;
	// cout << "_BC_RHS = " << endl;
	// cout << _BC_RHS << endl;
	// cout << "-------------------------------" << endl;
}

double FiniteVolume::max(double a, double b)
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


#define _FINITEVOLUME_CPP
#endif
