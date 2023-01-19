#include "md.h"

#ifndef _DYNAMICS_PARAM_H
#define _DYNAMICS_PARAM_H

class DynamicsParam {
  int nstep,log_freq,tra_freq;
  double dt,tau_t,ref_t;

  void clean(void);
  void copy(const DynamicsParam& a);

public:
  DynamicsParam(void);
  ~DynamicsParam(void);
  DynamicsParam(const DynamicsParam& a);
  DynamicsParam& operator = (const DynamicsParam& a);

  void     setup(void);
  int      get_nstep(void);
  int      get_log_freq(void);
  int      get_tra_freq(void);
  double   get_dt(void);
  double   get_ref_t(void);
  double   get_tau_t(void);
  void     set_nstep(int nstep);
  void     set_log_freq(int log_freq);
  void     set_tra_freq(int tra_freq);
  void     set_dt(double dt);
  void     set_ref_t(double ref_t);
  void     set_tau_t(double tau_t);
};
#endif
