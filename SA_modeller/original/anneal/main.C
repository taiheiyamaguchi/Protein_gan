#include "md.h"

int main(int ac,char **av)
{
  Molecule a;
  Potential pot;
  DynamicsParam dp;
  int i,natom,ndim=3,npair,**list;
  double *mass;
  FILE *fin,*flog,*fpdb;

// Usage
  if(ac != 8) {
    fprintf(stderr,"%s nsteps output.log log_freq noe.tbl tra.pdb tra_freq final.pdb\n",*av);
    exit(1);
  }

// Read Distance restraints
  fin=fopen(*(av+4),"r");
  pot.read(fin);
  npair=pot.get_npair();
  fclose(fin);
  printf("Number of restraints: %d\n",npair);

// Find number of atoms
  list=pot.get_list();
  natom=0;
  for(i=0;i<npair;i++) {
    if(natom < list[i][0]) {
      natom=list[i][0];
    }
    if(natom < list[i][1]) {
      natom=list[i][1];
    }
  }
  natom++;
  printf("Number of atoms: %d\n",natom);
  printf("Number of dimensions: %d\n",ndim);

// Set up molecule
  a.setup(natom,ndim);
  mass=a.get_mass();
  for(i=0;i<natom;i++) {
    mass[i]=12.0;
  }
  a.random_x(20.0);
  a.random_v(300.0);

  dp.set_nstep(atoi(*(av+1)));
  dp.set_log_freq(atoi(*(av+3)));
  dp.set_tra_freq(atoi(*(av+6)));
  dp.set_dt(0.001);
  dp.set_ref_t(300.0);
  dp.set_tau_t(0.01);

  flog=fopen(*(av+2),"w");
  fpdb=fopen(*(av+5),"w");
  a.dynamics(pot,dp,flog,fpdb);
  fclose(flog);
  fclose(fpdb);

  fpdb=fopen(*(av+7),"w");
  a.writePDB(fpdb);
  fclose(fpdb);
  
  return 0;
}
