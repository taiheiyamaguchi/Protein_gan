#include "md.h"

Potential::Potential(void)
{
  npair=0;
  list=NULL;
  dist=NULL;
  fc=NULL;
  epot=0.0;
}

Potential::~Potential(void)
{
  this->clean();
}

Potential::Potential(const Potential& pot)
{
  this->copy(pot);
}

Potential& Potential::operator = (const Potential& pot)
{
  if(this == &pot)
    return *this;

  this->clean();
  this->copy(pot);
  return *this;
}

void Potential::clean(void)
{
  FREE_MATRIX(list,npair);
  FREE_MATRIX(dist,npair);
  FREE_ARRAY(fc);
}

void Potential::copy(const Potential& pot)
{
  int i,j;

  this->setup(pot.npair);

  for(i=0;i<npair;i++) {
    list[i][0]=pot.list[i][0];
    list[i][1]=pot.list[i][1];
    dist[i][0]=pot.dist[i][0];
    dist[i][1]=pot.dist[i][1];
    dist[i][2]=pot.dist[i][2];
    fc[i]=pot.fc[i];
  }

  epot=pot.epot;
}

void Potential::setup(int npair1)
{
  npair=npair1;

  ALLOCATE_INT_MATRIX(list,npair,2);
  ALLOCATE_DOUBLE_MATRIX(dist,npair,3);
  ALLOCATE_DOUBLE_ARRAY(fc,npair);
}

void Potential::read(FILE *fp)
{
  int i0,i1,num,i;
  float d0,d1,d2,d3;

  num=0;
  while(1) {
    if(6 != fscanf(fp,"%d %d %f %f %f %f",&i0,&i1,&d0,&d1,&d2,&d3)) {
      break;
    } else {
      num++;
    }
  }
  this->setup(num);
  rewind(fp);
  for(i=0;i<npair;i++) {
    fscanf(fp,"%d %d %f %f %f %f",&i0,&i1,&d0,&d1,&d2,&d3);
    list[i][0]=i0;
    list[i][1]=i1;
    dist[i][0]=(double)d0;
    dist[i][1]=(double)d1;
    dist[i][2]=(double)d2;
    fc[i]=(double)d3;
  }
  epot=0.0;
}

void Potential::calc_force(double **x,double **f,int natom,int ndim)
{
  int i,j,at0,at1;
  double *d,r0,r1,delta,fac,r,e,de;

  for(i=0;i<natom;i++) {
    for(j=0;j<ndim;j++) {
      f[i][j]=0.0;
    }
  }

  ALLOCATE_DOUBLE_ARRAY(d,ndim);
  epot=0.0;
  for(i=0;i<npair;i++) {
    at0=list[i][0];
    at1=list[i][1];
    r0=dist[i][0];
    r1=dist[i][1];
    delta=dist[i][2];
    fac=dist[i][3];
    r=0.0;
    for(j=0;j<ndim;j++) {
      d[j]=x[at0][j]-x[at1][j];
      r+=SQR(d[j]);
    }
    r=sqrt(r);
    if(r < r0-delta) {
      e=-fac*delta*(r-r0+0.5*delta);
      de=-fac*delta;
    } else if(r < r0) {
      e=0.5*fac*SQR(r-r0);
      de=fac*(r-r0);
    } else if(r < r1) {
      e=0.0;
      de=0.0;
    } else if(r < r1+delta) {
      e=0.5*fac*SQR(r-r1);
      de=fac*(r-r1);
    } else {
      e=fac*delta*(r-r1-0.5*delta);
      de=fac*delta;
    }
    for(j=0;j<ndim;j++) {
      f[at0][j]-=de*d[j]/r;
      f[at1][j]+=de*d[j]/r;
    }
    epot+=e;
  }
  FREE_ARRAY(d);
}

int Potential::get_npair(void)
{
  return npair;
}

int** Potential::get_list(void)
{
  return list;
}

double** Potential::get_dist(void)
{
  return dist;
}

double* Potential::get_fc(void)
{
  return fc;
}

double Potential::get_epot(void)
{
  return epot;
}
