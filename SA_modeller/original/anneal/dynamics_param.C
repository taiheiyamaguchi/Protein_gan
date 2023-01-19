#include "md.h"

DynamicsParam::DynamicsParam(void)
{
  nstep=0;
  log_freq=0;
  tra_freq=0;
  dt=0.0;
  ref_t=0.0;
  tau_t=0.0;
}

DynamicsParam::~DynamicsParam(void)
{
  this->clean();
}

DynamicsParam::DynamicsParam(const DynamicsParam& dp)
{
  this->copy(dp);
}

DynamicsParam& DynamicsParam::operator = (const DynamicsParam& dp)
{
  if(this == &dp)
    return *this;

  this->clean();
  this->copy(dp);
  return *this;
}

void DynamicsParam::clean(void)
{
}

void DynamicsParam::copy(const DynamicsParam& dp)
{
  this->setup();

  nstep=dp.nstep;
  log_freq=dp.log_freq;
  tra_freq=dp.tra_freq;
  dt=dp.dt;
  ref_t=dp.ref_t;
  tau_t=dp.tau_t;
}

void DynamicsParam::setup(void)
{
}

int DynamicsParam::get_nstep(void)
{
  return nstep;
}

int DynamicsParam::get_log_freq(void)
{
  return log_freq;
}

int DynamicsParam::get_tra_freq(void)
{
  return tra_freq;
}

double DynamicsParam::get_dt(void)
{
  return dt;
}

double DynamicsParam::get_ref_t(void)
{
  return ref_t;
}

double DynamicsParam::get_tau_t(void)
{
  return tau_t;
}

void DynamicsParam::set_nstep(int nstep1)
{
  nstep=nstep1;
}

void DynamicsParam::set_log_freq(int log_freq1)
{
  log_freq=log_freq1;
}

void DynamicsParam::set_tra_freq(int tra_freq1)
{
  tra_freq=tra_freq1;
}

void DynamicsParam::set_dt(double dt1)
{
  dt=dt1;
}

void DynamicsParam::set_ref_t(double ref_t1)
{
  ref_t=ref_t1;
}

void DynamicsParam::set_tau_t(double tau_t1)
{
  tau_t=tau_t1;
}
