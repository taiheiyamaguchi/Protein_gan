#include "md.h"

Molecule::Molecule(void)
{
  natom=0;
  ndim=0;
  x=NULL;
  v=NULL;
  f=NULL;
  mass=NULL;
  ekin=0.0;
}

Molecule::Molecule(int natom)
{
  this->setup(natom,3);
}

Molecule::Molecule(int natom,int ndim)
{
  this->setup(natom,ndim);
}

Molecule::~Molecule(void)
{
  this->clean();
}

Molecule::Molecule(const Molecule& a)
{
  this->copy(a);
}

Molecule& Molecule::operator = (const Molecule& a)
{
  if(this == &a)
    return *this;

  this->clean();
  this->copy(a);
  return *this;
}

void Molecule::clean(void)
{
  FREE_MATRIX(x,natom);
  FREE_MATRIX(v,natom);
  FREE_MATRIX(f,natom);
  FREE_ARRAY(mass);
}

void Molecule::copy(const Molecule& a)
{
  int i,j;

  this->setup(a.natom,a.ndim);

  for(i=0;i<natom;i++) {
    mass[i]=a.mass[i];
    for(j=0;j<ndim;j++) {
      x[i][j]=a.x[i][j];
      v[i][j]=a.v[i][j];
      f[i][j]=a.f[i][j];
    }
  }

  ekin=a.ekin;
}

void Molecule::setup(int natom1,int ndim1)
{
  natom=natom1;
  ndim=ndim1;

  ALLOCATE_DOUBLE_MATRIX(x,natom,ndim);
  ALLOCATE_DOUBLE_MATRIX(v,natom,ndim);
  ALLOCATE_DOUBLE_MATRIX(f,natom,ndim);
  ALLOCATE_DOUBLE_ARRAY(mass,natom);
}

void Molecule::random_x(double size)
{
  int i,j;

  srand48(time(NULL));	/* Use epoch time as random seed */
  for(i=0;i<natom;i++) {
    for(j=0;j<ndim;j++) {
      x[i][j]=drand48()*size;
    }
  }
}

void Molecule::random_v(double temp)
{
  int i,j;
  double *r,x1,x2,y1,y2,fac;

  srand48(time(NULL));	/* Use epoch time as random seed */

  ALLOCATE_DOUBLE_ARRAY(r,natom*ndim+1);
  for(i=0;i<natom*ndim;i+=2) {
    x1=drand48();
    x2=drand48();
    y1=sqrt(-2.0*log(x1))*cos(2.0*PI*x2);
    y2=sqrt(-2.0*log(x1))*sin(2.0*PI*x2);
    r[i  ]=y1;
    r[i+1]=y2;
  }

  for(i=0;i<natom;i++) {
    fac=sqrt(BOLTZ*temp/mass[i]);
    for(j=0;j<ndim;j++) {
      v[i][j]=r[i*ndim+j]*fac;
    }
  }
  FREE_ARRAY(r);
}

void Molecule::calc_ekin(void)
{
  int i,j;

  ekin=0.0;
  for(i=0;i<natom;i++) { 
    for(j=0;j<ndim;j++) { 
      ekin+=mass[i]*SQR(v[i][j]);
    }
  }
  ekin*=0.5;
}

void Molecule::dynamics(Potential pot,DynamicsParam dparam,
                        FILE *flog,FILE *ftra)
{
  int nstep,i,j,k,log_freq,tra_freq,nmodel=0;
  double dt,dti,epot,temp,tau_t,ref_t,sca;

  nstep=dparam.get_nstep();
  dt=dparam.get_dt();
  dti=dt*20.455;
  log_freq=dparam.get_log_freq();
  tra_freq=dparam.get_tra_freq();
  tau_t=dparam.get_tau_t();
  ref_t=dparam.get_ref_t();
//
// Velocity Verlet
//
  pot.calc_force(x,f,natom,ndim);
  for(k=0;k<nstep;k++) {
    for(i=0;i<natom;i++) {
      for(j=0;j<ndim;j++) {
        v[i][j]+=0.5*f[i][j]/mass[i]*dti;
        x[i][j]+=v[i][j]*dti;
      }
    }
    pot.calc_force(x,f,natom,ndim);
    for(i=0;i<natom;i++) {
      for(j=0;j<ndim;j++) {
        v[i][j]+=0.5*f[i][j]/mass[i]*dti;
      }
    }
//
// Berendsen thermostat
//
    if(ref_t > 0) {
      this->calc_ekin();
      temp=ekin*2.0/((double)(ndim*natom)*BOLTZ);
      sca=sqrt(1.0+dt/tau_t*(ref_t/temp-1.0));
      for(i=0;i<natom;i++) {
        for(j=0;j<ndim;j++) {
          v[i][j]*=sca;
        }
      }
    }
//
// Print log
//
    if((k+1) % log_freq == 0) {
      this->calc_ekin();
      temp=ekin*2.0/((double)(ndim*natom)*BOLTZ);
      epot=pot.get_epot();
      fprintf(flog,"Step %8d TIME %.3f EPOT %e EKIN %e ETOT %e TEMP %.4f\n",
        k+1,dt*(k+1),epot,ekin,epot+ekin,temp);
    }
//
// Save trajectory
//
    if((k+1) % tra_freq == 0) {
      nmodel++;
      fprintf(ftra,"MODEL %8d\n",nmodel);
      this->writePDB(ftra);
      fprintf(ftra,"ENDMDL\n");
    }
  }
}

void Molecule::writePDB(FILE *fp)
{
  int i;

  for(i=0;i<natom;i++) {
    fprintf(fp,"%-6s%5d %-4s %-4s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
      "ATOM",i+1,"CA","GLY",i+1,x[i][0],x[i][1],x[i][2],1.0,0.0);
//    "HETATM",i+1,"K+","K+",i+1,x[i][0],x[i][1],x[i][2],1.0,0.0);
  }
}

double** Molecule::get_x(void)
{
  return x;
}

double** Molecule::get_v(void)
{
  return v;
}

double** Molecule::get_f(void)
{
  return f;
}

double* Molecule::get_mass(void)
{
  return mass;
}

double Molecule::get_ekin(void)
{
  return ekin;
}
