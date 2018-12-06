#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "tools.h"

int startswith(const char *pre, const char *str)
{
  return strncmp(pre, str, strlen(str)) == 0;
}

/* merge multiple space into one space and remove endline whitespace */
void strip_str(char *str)
{
  int i=0, j;
  int len;
  len = strlen(str);
  while(i<len) {
    if(str[i]==' ' && (str[i+1]==' ' || str[i-1]==' ')) {
      for(j=i; j<len; j++)
        str[j]=str[j+1];
      len--;
    }
    else if(str[i] == ' ' && str[i+1] == '\n')
      len--;
    else
      i++;
  }
}

//get io file directory
char *get_dir(char *fnm)
{
  int i;
  char* curdir;
  size_t size = strlen(fnm);

  for (i = (size - 1); (i >= 0) && (fnm[i] != '/'); --i);
  if (fnm[i] == '/') { 
    curdir = (char*) malloc((i + 1) * sizeof(char));
    if (curdir != NULL) {
      strncpy(curdir, fnm, (size_t) i);
      curdir[i] = '\0';
      return curdir;
    }
    else 
      fprintf(stderr,"Memory allocation error.");
  }
  else
    fprintf(stderr,"There was no DIR_SEPARATOR in LOC.");
  
  curdir = ".";
  return curdir;
}

//vector, matrix operations
void init_matrix(int nrow, int ncol, float **(*pmat))
{
  int i, j;

  *pmat = (float **) malloc (nrow * sizeof(float *));
  for(i=0; i<nrow; i++)
    (*pmat)[i] = (float *) malloc (ncol * sizeof(float));

  for(i=0; i<nrow; i++)
    for(j=0; j<ncol; j++)
      (*pmat)[i][j] = 0.;
}

void init_state(t_state *md_state, int nr)
{
  md_state->x = malloc (nr * sizeof(fvec));
  md_state->v = malloc (nr * sizeof(fvec));
  md_state->f = malloc (nr * sizeof(fvec));
}

void matrix_scale(int nrow, int ncol, float **ma, float sFactor)
{
  int i, j;

  for(i=0; i<nrow; i++)
    for(j=0; j<ncol; j++)
      ma[i][j] *= sFactor;
}

void matrix_shift(int nrow, int ncol, float **ma, float sFactor)
{
  int i, j;

  for(i=0; i<nrow; i++)
    for(j=0; j<ncol; j++)
      ma[i][j] -= sFactor;
}

void matrix_sum(int nrow, int ncol, float **ma1, float **ma2)
{
  int i, j;

  for(i=0; i<nrow; i++)
    for(j=0; j<ncol; j++)
      ma1[i][j] = ma1[i][j] + ma2[i][j];
}

float vec_dot(fvec v1, fvec v2)
{
  return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}


//file operation
void parse_comm_args(int *argc, char *argv[], t_mdPara *md_para,
                    char *fn_in, char *fn_out)
{
  int nCubic, nPerEdge;
  char *fn, *tmpstring = NULL;
  size_t tmplen = 0;
  FILE *fp;
  char lineEle1[50], lineEle2;

  if(*argc < 2) {
    fprintf(stderr, "no enough arguments\n");
    exit(1);
  }
  else if(*argc > 2) {
    fprintf(stderr, "Warning: Arguments after %s are to be omitted\n", argv[1]);
  }

  fn = argv[1];
  fp = fopen(fn, "r");
  if(fp == NULL) {
    fprintf(stderr, "md parameter file does not exist\n");
    exit(1);
  }

  //default value
  md_para->integrator = 1;
  md_para->dt = 0.005;
  md_para->mass = 1;
  md_para->density = 0.8;
  md_para->epsilon = 1;
  md_para->sigma = 1;
  md_para->rCutoff = 2.5;
  //parse the input parameter file
  while (!feof(fp)) {
    getline(&tmpstring, &tmplen, fp);
    strip_str(tmpstring);
    if(startswith(tmpstring, "#") || startswith(tmpstring, "\n"))
      continue;
    else if(startswith(tmpstring, "input"))
      sscanf(tmpstring, "%s %s %s", lineEle1, &(lineEle2), fn_in);
    else if(startswith(tmpstring, "output"))
      sscanf(tmpstring, "%s %s %s", lineEle1, &(lineEle2), fn_out);
    else if(startswith(tmpstring, "integrator"))
      sscanf(tmpstring, "%s %s %d", lineEle1, &(lineEle2), &(md_para->integrator));
    else if(startswith(tmpstring, "dt"))
      sscanf(tmpstring, "%s %s %f", lineEle1, &(lineEle2), &(md_para->dt));
    else if(startswith(tmpstring, "nSteps"))
      sscanf(tmpstring, "%s %s %d", lineEle1, &(lineEle2), &(md_para->nSteps));
    else if(startswith(tmpstring, "mass"))
      sscanf(tmpstring, "%s %s %f", lineEle1, &(lineEle2), &(md_para->mass));
    else if(startswith(tmpstring, "density"))
      sscanf(tmpstring, "%s %s %f", lineEle1, &(lineEle2), &(md_para->density));
    else if(startswith(tmpstring, "natoms"))
      sscanf(tmpstring, "%s %s %d", lineEle1, &(lineEle2), &(md_para->natoms));
    else if(startswith(tmpstring, "temperature"))
      sscanf(tmpstring, "%s %s %f", lineEle1, &(lineEle2), &(md_para->temperature));
    else if(startswith(tmpstring, "epsilon"))
      sscanf(tmpstring, "%s %s %f", lineEle1, &(lineEle2), &(md_para->epsilon));
    else if(startswith(tmpstring, "sigma"))
      sscanf(tmpstring, "%s %s %f", lineEle1, &(lineEle2), &(md_para->sigma));
    else if(startswith(tmpstring, "rCutoff"))
      sscanf(tmpstring, "%s %s %f", lineEle1, &(lineEle2), &(md_para->rCutoff));
    else {
      fprintf(stderr, "parameter %s is not defined\n", lineEle1);
      exit(0);
    }
  }

  md_para->vCutoff = 4 *  md_para->epsilon *
                    (pow((md_para->sigma/md_para->rCutoff), 12)
                   - pow((md_para->sigma/md_para->rCutoff), 6));
  nCubic = (md_para->natoms)/4;
  nPerEdge = (int) (pow(nCubic, 1./3) + 0.5);
  md_para->boxLength = pow((4. * md_para->mass/ md_para->density), 1./3) * nPerEdge;
  md_para->nPerEdge = nPerEdge;
  md_para->natoms = 4 * pow(nPerEdge, 3); //when nPerEdge is not a integrator;

  fclose(fp);
}

void write_xyz(FILE *fp, t_mdPara *md_para, float **x)
{
  int i, nr;

  nr = md_para->natoms;
  fprintf(fp, "%d\n", nr);
  fprintf(fp, "fcc lattice with density %f\n", md_para->density);
  for(i=0; i<nr; i++) {
    fprintf(fp, "H %lf %lf %lf\n", x[i][0], x[i][1], x[i][2]);
  }
  fprintf(fp, "\n");
}

void write_all(FILE *fp, t_mdPara *md_para, t_state *md_state)
{
  int i, nr;

  nr = md_para->natoms;
  fprintf(fp, "%d\n", nr);
  fprintf(fp, "kinE %lf potE %lf\n", md_state->kinE, md_state->potE);
  for(i=0; i<nr; i++) {
    fprintf(fp, "H %lf %lf %lf %lf %lf %lf\n",
              md_state->x[i][0], md_state->x[i][1], md_state->x[i][2],
              md_state->v[i][0], md_state->v[i][1], md_state->v[i][2]);
  }
}

void read_xyz(char *fn, t_state *md_state, t_mdPara *md_para)
{
  FILE *fp;
  char tmpchar, *tmpstring = NULL;
  size_t tmplen = 0;
  int i, nr;

  fp = fopen(fn, "r");
  if(fp == NULL) {
    fprintf(stderr, "input xyz file does not exist\n");
    exit(1);
  }
  getline(&tmpstring, &tmplen, fp);
  strip_str(tmpstring);
  sscanf(tmpstring, "%d", &nr);
  if(nr != md_para->natoms) {
    fprintf(stderr, "atom number in xyz file is not consistent with parameter file\n");
    exit(1);
  }

  getline(&tmpstring, &tmplen, fp);
  for(i=0; i<nr; i++) {
    if(!(feof(fp))) {
      getline(&tmpstring, &tmplen, fp);
      strip_str(tmpstring);
      sscanf(tmpstring, "%c %f %f %f", &tmpchar, &(md_state->x[i][0]),
             &(md_state->x[i][1]), &(md_state->x[i][2]));
    }
  }
  if(i != nr) {
    fprintf(stderr, "number of coordinates is not equal to atom numbers\n");
    exit(1);
  }

  fclose(fp);
}

int read_traj(FILE *fp, t_state *md_state, t_mdPara *md_para)
{
  char tmpchar, tmpstring[4], *tmpline = NULL;
  size_t tmplen = 0;
  int i, nr;

  if(feof(fp)) {
    fprintf(stderr, "finishing reading trajectory file\n");
    return 0;
  }

  nr = md_para->natoms;
  getline(&tmpline, &tmplen, fp);
  getline(&tmpline, &tmplen, fp);
  strip_str(tmpline);
  sscanf(tmpline, "%s %f %s %f", tmpstring, &(md_state->kinE), tmpstring,
		&(md_state->potE));
  for(i=0; i<nr; i++) {
    getline(&tmpline, &tmplen, fp);
    strip_str(tmpline);
    sscanf(tmpline, "%c %f %f %f %f %f %f", &tmpchar,
           &(md_state->x[i][0]), &(md_state->x[i][1]), &(md_state->x[i][2]),
           &(md_state->v[i][0]), &(md_state->v[i][1]), &(md_state->v[i][2]));
  }
  if(i != nr) {
    fprintf(stderr, "number of coordinates is not equal to atom numbers\n");
    exit(1);
  }
  return 0;
}

//related to simulation
float calc_r2(fvec x)
{
  return (pow(x[0], 2) + pow(x[1], 2) + pow(x[2], 2));
}

//generating random numbers from normal distribution
float randn(float mu, float sigma)
{
  float U1, U2, W, mult;
  static float X1, X2;
  static int call = 0;

  if(call == 1) {
    call = !call;
    return (mu + sigma * (float) X2);
  }
  do {
    U1 = -1 + ((float) rand () / RAND_MAX) * 2;
    U2 = -1 + ((float) rand () / RAND_MAX) * 2;
    W = pow(U1, 2) + pow(U2, 2);
  } while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
  call = !call;

  return (mu + sigma * (float) X1);
}

void gen_velocity(t_mdPara *md_para, t_state *md_state)
{
  int i, j, nr;
  float temperature, mass, epsilon, sigma;

  srand((int)time(NULL));
  nr = md_para->natoms;
  mass = md_para->mass;
  epsilon = md_para->epsilon;
  temperature = md_para->temperature;
  sigma = mass * epsilon / temperature;

  for(i=0; i<nr; i++)
    for(j=0; j<3; j++)
      md_state->v[i][j] = randn(0, sigma);
  calc_kinE(md_para, md_state);
}

void calc_kinE(t_mdPara *md_para, t_state *md_state)
{
  int i, j;
  float mass, kinE;

  mass = md_para->mass;
  kinE = 0;
  for(i=0; i<md_para->natoms; i++)
    for(j=0; j<3; j++)
      kinE += 0.5 * mass * (md_state->v)[i][j] * (md_state->v)[i][j];
  md_state->kinE = kinE;
}

void get_vec_pbc(float length, fvec x1, fvec x2, fvec vecij)
{
  int i;

  for(i=0; i<3; i++) {
    vecij[i] = x2[i] - x1[i];
    vecij[i] = vecij[i] - length * round(vecij[i]/length);
  }
}

float calc_lj(float c12, float c6, float vCut, int indi, int indj,
              fvec vec_ij, float r_sq, fvec *f)
{
  int i;
  float r6, r12, potE, f1, f2;

  r6 = pow(r_sq, 3);
  r12 = pow(r6, 2);
  f1 = 12 * c12 / (r12 * r_sq);
  f2 = 6 * c6 / (r6 * r_sq);

  for(i=0; i<3; i++) {
    f[indi][i] -= (f1 - f2) * vec_ij[i];
    f[indj][i] += (f1 - f2) * vec_ij[i];
  }
  potE = c12 / r12 - c6 / r6 - vCut;

  return potE;
}

float calc_force(t_mdPara *md_para, fvec *x, fvec *f)
{
  int i, j, k, nr;
  float mass, C12, C6, rij_sq, rCut_sq, length, potE;
  fvec vec_ij;

  nr = md_para->natoms;
  length = md_para->boxLength;
  mass = md_para->mass;
  rCut_sq = pow(md_para->rCutoff, 2);
  C12 = 4 * md_para->epsilon * pow(md_para->sigma, 12);
  C6 = 4 * md_para->epsilon * pow(md_para->sigma, 6);

  potE = 0;
  for(i=0; i<nr; i++)
    for(k=0; k<3; k++)
      f[i][k] = 0;

  for(i=0; i<nr-1; i++) {
    for(j=i+1; j<nr; j++) {
      get_vec_pbc(length, x[i], x[j], vec_ij);
      rij_sq = vec_ij[0]*vec_ij[0] + vec_ij[1]*vec_ij[1] + vec_ij[2]*vec_ij[2];
      if(rij_sq < rCut_sq) {
        if(rij_sq < 0.56) {
          fprintf(stderr, "atom %d and %d has distance %f < 0.56, overlap\n", i, j,rij_sq);
          exit(1);
        }
	potE += calc_lj(C12, C6, md_para->vCutoff, i, j, vec_ij, rij_sq, f);
      }
    }
  }
  return potE;
}

void update(t_mdPara *md_para, t_state *md_state)
{
  int i, j, nr;
  float dt, dt2, mass, length, length2;
  fvec *f_old, *f_new;

  nr = md_para->natoms;
  dt = md_para->dt;
  dt2 = dt * dt;
  mass = md_para->mass;
  length = md_para->boxLength;
  length2 = length/2;
  f_old = md_state->f;
  f_new = malloc(nr*sizeof(fvec));

  for(i=0; i<nr; i++) {
    for(j=0; j<3; j++) {
      md_state->x[i][j] += dt * md_state->v[i][j] + 0.5 * dt2 * f_old[i][j]/mass;
      if(md_state->x[i][j] > length2)
        md_state->x[i][j] = md_state->x[i][j] - length;
      else if (md_state->x[i][j] < -1*length2)
        md_state->x[i][j] = md_state->x[i][j] + length;
    }
  }
  md_state->potE = calc_force(md_para, md_state->x, f_new);
  for(i=0; i<nr; i++) {
    for(j=0; j<3; j++) {
      md_state->v[i][j] += 0.5 * dt * (f_old[i][j] + f_new[i][j]) / mass;
      md_state->f[i][j] = f_new[i][j];
    }
  }
  calc_kinE(md_para, md_state);

  free(f_new);
}

