#include "md.h"

#ifndef _MOLECULE_H
#define _MOLECULE_H

class Molecule {
  int natom;
  int ndim;
  double **x;
  double **v;
  double **f;
  double *mass;
  double ekin;

  void clean(void);
  void copy(const Molecule& a);

public:
  Molecule(void);
  Molecule(int natom);
  Molecule(int natom,int ndim);
  ~Molecule(void);
  Molecule(const Molecule& a);
  Molecule& operator = (const Molecule& a);

  void     setup(int natom,int ndim);
  void     random_x(double size);
  void     random_v(double temp);
  void     calc_ekin(void);
  void     dynamics(Potential pot,DynamicsParam dparam,FILE *flog,FILE *fpdb);
  void     writePDB(FILE *fp);
  double** get_x(void);
  double** get_v(void);
  double** get_f(void);
  double*  get_mass(void);
  double   get_ekin(void);
};
#endif
