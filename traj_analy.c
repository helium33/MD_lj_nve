//main function for md simulation
//Shanshan Wu
//12/03/2018
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tools.h"
#define Nbins 170
#define dr 0.02

//radius distribution function
int get_index(float m)
{
  return((int) (m * 10000) / 10000 +1);
}

float get_v(float r)
{
  return (4*M_PI*dr*r*r);
}

void gr_1frame(t_state *md_state, t_mdPara *md_para, int *gr)
{
  int i, j, nr, grIndex;
  float length, rij;
  fvec vec_ij;

  nr = md_para->natoms;
  length = md_para->boxLength;
  for(i=0; i<nr-1; i++) {
    for(j=i+1; j<nr; j++) {
      get_vec_pbc(length, md_state->x[i], md_state->x[j], vec_ij);
      rij = sqrt(calc_r2(vec_ij));
      grIndex = get_index(rij/0.02);
      if(grIndex < Nbins)
        gr[grIndex]++;
    }
  }
}

//translational order parameter
float transOP_1frame(t_state *md_state, t_mdPara *md_para, fvec vec_k, fvec vec0)
{
  int i, j;
  fvec vec;
  float trans=0;

  for(i=0; i<md_para->natoms; i++) {
    for(j=0; j<3; j++)
      vec[j] = md_state->x[i][j] - vec0[j];
    trans += cos(vec_dot(vec_k, vec));
  }

  return trans;
}

int main(int argc, char *argv[])
{
  t_mdPara md_para;
  t_state md_state;
  char *dir, fn_in[50], fn_out[50], fn_ene[50], fn_gr[50], fn_trans[50];
  FILE *fp_traj, *fp_ene, *fp_gr, *fp_trans;
  int i, istep, natoms;
  int gr[Nbins]; //r = 0-3.4, 0.02 bin size, half box size
  float volume, gr_coef, cellLength, transOP;
  fvec vec_k, vec0; //k vector for translational OP

  parse_comm_args(&argc, argv, &md_para, fn_in, fn_out);
  dir = get_dir(fn_in);
  natoms = md_para.natoms;
  cellLength = (md_para.boxLength) / (md_para.nPerEdge);
//Initialization and setup
  init_state(&md_state, md_para.natoms); 
  sprintf(fn_ene, "%s/traj_energy.dat", dir);
  sprintf(fn_gr, "%s/traj_gr.dat", dir);
  sprintf(fn_trans, "%s/traj_trans.dat", dir);
  fp_traj = fopen(fn_out, "r");
  if(fp_traj == NULL) {
    fprintf(stderr, "trajectory file does not exist\n");
    exit(1);
  }
  fp_ene = fopen(fn_ene, "w");
  fp_gr = fopen(fn_gr, "w");
  fp_trans = fopen(fn_trans, "w");
  if(fp_ene == NULL || fp_gr == NULL || fp_trans == NULL) {
    fprintf(stderr, "Can not open file to write\n");
    exit(1);
  }
  fprintf(stderr, "Doing trajectory analysis\n");
  fprintf(stderr, "1. total energy vs time and kinetic energy vs time\n");
  fprintf(stderr, "2. radial distribution function \n");
  fprintf(stderr, "3. translational order parameter\n");

  for(i=0; i<Nbins; i++)
    gr[i] = 0;
  
  read_traj(fp_traj, &md_state, &md_para);
  for(i=0; i<3; i++) {
    vec_k[i] = -2 * M_PI / cellLength;
    vec0[i] = md_state.x[0][i];
  }
  for(istep=0; istep<md_para.nSteps-1; istep++) {
// output energy 
    fprintf(fp_ene, "%d %f %f %f \n", istep, md_state.kinE, md_state.potE, 
            md_state.kinE+md_state.potE);
//radius distribution
    if(istep > 100)
      gr_1frame(&md_state, &md_para, gr);
//translational order parameter
    transOP = transOP_1frame(&md_state, &md_para, vec_k, vec0);
    fprintf(fp_trans, "%d %f \n", istep, transOP/natoms);
    read_traj(fp_traj, &md_state, &md_para);
  }
  gr_coef = 1/(md_para.density * natoms / 2 * (istep-100));
  for(i=1; i<Nbins; i++) {
    volume = get_v(i*dr);
    fprintf(fp_gr, "%f %f 1\n", i*dr-dr, gr[i] * gr_coef / volume);
  }

  fclose(fp_ene);
  fclose(fp_gr);
  fclose(fp_trans);
  fclose(fp_traj);

  return 0;
}

