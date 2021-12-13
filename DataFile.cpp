#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include <fstream>
#include <iostream>
#include "toml.hpp"
#include <math.h>

// Constructeur
DataFile::DataFile(std::string file_name)
{
  // Lecture du fichier de données
  auto config = toml::parse(file_name);

  // Other
  const auto& other = toml::find(config, "other");
  _results = toml::find<std::string>(other, "results");
  _lambdapv = toml::find<double>(other, "lambdapv");
  _ppv = toml::find<double>(other, "ppv");
  _cpv = toml::find<double>(other, "cpv");
  _Aref = toml::find<double>(other, "Aref");
  _Ta =  toml::find<double>(other, "Ta");
  _Lm = toml::find<double>(other, "Lm");
  _rhov =toml::find<double>(other, "rhov");
  _rhop = toml::find<double>(other, "rhop");
  _T0 = toml::find<double>(other, "T0");
  
  _sigma = _lambdapv/(_ppv*_cpv);
  system(("rm -r ./" + _results).c_str());
  system(("mkdir -p ./" + _results).c_str());

  // Time
  const auto& time = toml::find(config, "time");
  _t0 = toml::find<double>(time, "t0");
  _tfinal = toml::find<double>(time, "tfinal");
  _dt = toml::find<double>(time, "dt");
  _scheme = toml::find<std::string>(time, "scheme");

  // Space
  const auto& space = toml::find(config, "space");
  _xmin = toml::find<double>(space, "xmin");
  _xmax = toml::find<double>(space, "xmax");
  _dx = toml::find<double>(space, "dx");
  _ymin = toml::find<double>(space, "ymin");
  _ymax = toml::find<double>(space, "ymax");
  _dy = toml::find<double>(space, "dy");

  // Boundary conditions
  // const auto& BC = toml::find(config, "BC");
  // _LBC = toml::find<std::string>(BC, "LeftBoundCond");
  // _RBC = toml::find<std::string>(BC, "RightBoundCond");
  // _DBC = toml::find<std::string>(BC, "DownBoundCond");
  // _UBC = toml::find<std::string>(BC, "UpBoundCond");

  std::cout << "--------------------------------------------------" << std::endl;
  std::cout << "-------------- Adapt dt, dx and dy ---------------" << std::endl;
  std::cout << "-------------- xmax = xmin + Nx*dx -----------" << std::endl;
  std::cout << "-------------- ymax = ymin + Ny*dy -----------" << std::endl;
  std::cout << "-------------- tfinal = t0 + nb_it*dt ------------" << std::endl;
  // Calcul de _Nx et adaptation de _dx pour que (xmax - xmin) = (Nx+1)*dx
  _Nx = int(ceil((_xmax-_xmin)/_dx));
  _dx = (_xmax-_xmin)/(_Nx);
  // Calcul de _Ny et adaptation de _dy pour que (ymax - ymin) = (Ny+1)*dy
  _Ny = int(ceil((_ymax-_ymin)/_dy));
  _dy = (_ymax-_ymin)/(_Ny);

  // Calcul du pas de temps en fonction de la CFL
  // if (_scheme == "ExplicitEuler")
  // {
  //   if (_dt > 0.95*pow(_dy,2)*pow(_dx,2)/(2*_sigma*(pow(_dx,2)+pow(_dy,2))))
  //   {
  //     _dt=0.95*pow(_dy,2)*pow(_dx,2)/(2*_sigma*(pow(_dx,2)+pow(_dy,2)));
  //     std::cout << "The time step is fixed with the CFL condition: dt = " << _dt << "." << std::endl;
  //   }
  // }

  // Calcul du nombre d'itérations en temps
  // Adapter le pas de temps pour avoir _tfinal = _t0 + nb_iterations*_dt
  _nb_iterations = int(ceil((_tfinal-_t0)/_dt)); // Définition du nombre d’itérations(arrondi au plus grand)
  _dt = (_tfinal-_t0) / _nb_iterations;
  std::cout << "nb_iterations=" << _nb_iterations << " " << "_dt = " << _dt << " _dx = " << _dx << " _dy = " << _dy << " _Nx = " << _Nx << " _Ny = " << _Ny <<  std::endl;
  std::cout << "--------------------------------------------------" << std::endl;
}

#define _DATA_FILE_CPP
#endif
