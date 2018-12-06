//main function for md simulation
//Shanshan Wu                                                              
//12/03/2018                                                               
#include <stdio.h>                                                         
#include <stdlib.h>                                                        
#include <string.h>                                                        
#include <math.h>                                                          
#include "tools.h" 

int main(int argc, char * argv[])
{
  int iStep;
  t_mdPara md_para;   
  t_state md_state; 
  char fn_in[50], fn_out[50]; 
  FILE *fp_traj;

  //parse command line arguments                                             
  fprintf(stderr, "MD simulation of atomic liquid with L-J potential\n");
  parse_comm_args(&argc, argv, &md_para, fn_in, fn_out);
  fp_traj = fopen(fn_out, "w");
  init_state(&md_state, md_para.natoms);
  read_xyz(fn_in, &md_state, &md_para);
  gen_velocity(&md_para, &md_state);
  md_state.potE = calc_force(&md_para, md_state.x, md_state.f);
  write_all(fp_traj, &md_para, &md_state);
  for(iStep = 0; iStep < md_para.nSteps; iStep++) {
    update(&md_para, &md_state);
    write_all(fp_traj, &md_para, &md_state);
  }
 
  fclose(fp_traj);
  return 0;
}

