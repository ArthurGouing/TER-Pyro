#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>


class DataFile {
private:
  // Tous les paramètres contenus dans le fichier de données
  double _sigma;
  double _lambdapv;
  double _ppv;
  double _cpv;
  double _Aref, _Ta, _rhov, _rhop, _Lm, _T0;
  double _xmin;
  double _xmax;
  double _dx;
  int _Nx;
  double _ymin;
  double _ymax;
  double _dy;
  int _Ny;
  double _t0;
  double _tfinal;
  double _dt;
  int _nb_iterations;
  std::string _scheme;
  std::string _results;
  std::string _LBC, _RBC, _DBC, _UBC;
  // Nom du fichier de données
  const std::string _file_name;
  // Pour savoir si l'utilisateur a donné la valeur du paramètres
  // ou les paramètres par défaut doivent être utilisés
  bool _if_xmin;
  bool _if_xmax;
  bool _if_dx;
  bool _if_ymin;
  bool _if_ymax;
  bool _if_dy;
  bool _if_sigma;
  bool _if_t0;
  bool _if_tfinal;
  bool _if_dt;
  bool _if_scheme;
  bool _if_results;
  bool _if_LBC, _if_RBC, _if_DBC, _if_UBC;

public:
  // Constructeur
  DataFile(std::string file_name);

  const double & Get_xmin() const {return _xmin;};
  const double & Get_xmax() const {return _xmax;};
  const double & Get_dx() const {return _dx;};
  const int & Get_Nx() const {return _Nx;};
  const double & Get_ymin() const {return _ymin;};
  const double & Get_ymax() const {return _ymax;};
  const double & Get_dy() const {return _dy;};
  const int & Get_Ny() const {return _Ny;};
  const double & Get_ppv() const {return _ppv;};
  const double & Get_cpv() const {return _cpv;};
  const double & Get_Aref() const {return _Aref;};
  const double & Get_Ta() const {return _Ta;};
  const double & Get_rhov() const {return _rhov;};
  const double & Get_rhop() const {return _rhop;};
  const double & Get_Lm() const {return _Lm;};
  const double & Get_sigma() const {return _sigma;};
  const double & Get_t0() const {return _t0;};
  const double & Get_tfinal() const {return _tfinal;};
  const double & Get_dt() const {return _dt;};
  const double & Get_T0() const {return _T0;};
  const int & Get__nb_iterations() const {return _nb_iterations;};
  const std::string & Get_scheme() const {return _scheme;};
  const std::string & Get_results() const {return _results;};
  const std::string & Get_LBC() const {return _LBC;};
  const std::string & Get_RBC() const {return _RBC;};
  const std::string & Get_DBC() const {return _DBC;};
  const std::string & Get_UBC() const {return _UBC;};
};

#define _DATA_FILE_H
#endif
