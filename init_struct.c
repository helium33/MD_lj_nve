// Generate the initial configuration of atomic liquid for md simulation. 
//Shanshan Wu
//12/03/2018
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

void gen_fcc(t_mdPara *md_para, float **(*px))                             
{                                                                          
  int ix, iy, iz, i, j, a;                                                 
  int nr, nPerEdge;                                            
  float r_fcc[4][3] = {{0.25, 0.25, 0.25},                                 
                       {0.25, 0.75, 0.75},                                 
                       {0.75, 0.75, 0.25},                                 
                       {0.75, 0.25, 0.75}};
  fvec com = {0,0,0};

  nr = md_para->natoms;                                          
  nPerEdge = md_para->nPerEdge;;                              

  fprintf(stderr, "generate fcc lattice with %d atoms\n", nr);
  init_matrix(nr, 3, px); 

  i = 0;                                                                   
  for(ix=0; ix<nPerEdge; ix++) {                                           
    for(iy=0; iy<nPerEdge; iy++) {                                         
      for(iz=0; iz<nPerEdge; iz++) {                                       
        for(a=0; a<4; a++) {                                               
          (*px)[i][0] = r_fcc[a][0] + ix;                                  
          (*px)[i][1] = r_fcc[a][1] + iy;                                  
          (*px)[i][2] = r_fcc[a][2] + iz;                                  
          i++;                                                             
        }                                                                  
      }                                                                    
    }                                                                      
  }                                                                        
  matrix_scale(nr, 3, *px, md_para->boxLength/nPerEdge);
  for(j=0; j<3; j++) {
    for(i=0; i<nr; i++)
      com[j] += (*px)[i][j];
    com[j] /= nr;
  }
  matrix_shift(nr, 3, *px, com);                        
}    

int main(int argc, char* argv[])
{
  t_mdPara md_para; 
  float **x; 
  FILE *fp_xyz;
  char fn_in[50], fn_out[50];

//parse command line arguments
  parse_comm_args(&argc, argv, &md_para, fn_in, fn_out);
  
  fp_xyz = fopen(fn_in, "w");
  if(fp_xyz == NULL) {                                                         
    fprintf(stderr, "can not open file to write\n");                       
    exit(1);                                                               
  } 

  gen_fcc(&md_para, &x);
  write_xyz(fp_xyz, &md_para, x);

  fclose(fp_xyz);
  return 0;
}
