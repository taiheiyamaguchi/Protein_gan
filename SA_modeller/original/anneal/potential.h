#include "md.h"

#ifndef _POTENTIAL_H
#define _POTENTIAL_H

class Potential {
  int npair;
  int **list;		/* Atom pair list */
  double **dist;	/* dmin, dmax, delta */
  double *fc;		/* force constants in kcal/(mol A**2) */
  double epot;

  void clean(void);
  void copy(const Potential& pot);

public:
  Potential(void);
  ~Potential(void);
  Potential(const Potential& pot);
  Potential& operator = (const Potential& pot);

  void     setup(int npair);
  void     read(FILE *fp);
  void     calc_force(double **x,double **f,int natom,int ndim);
  int      get_npair(void);
  int**    get_list(void);
  double** get_dist(void);
  double*  get_fc(void);
  double   get_epot(void);
};
#endif
