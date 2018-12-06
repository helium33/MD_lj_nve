typedef float fvec[3];
 
typedef struct
{
  fvec *x;
  fvec *v;
  fvec *f;
  float potE;
  float kinE;
} t_state;

typedef struct 
{
  int integrator;
  int natoms;
  float mass;
  float density;
  float boxLength;
  int nPerEdge; //number of boxes per edge
  float epsilon;
  float sigma;
  float temperature;
  float rCutoff;
  float vCutoff;
  float dt;
  int nSteps;
} t_mdPara;

float calc_r2(fvec x);
float vec_dot(fvec x1, fvec x2);
char *get_dir(char *fnm);
void init_matrix(int nrow, int ncol, float **(*pmat));
void init_state(t_state *md_state, int nr);
void matrix_scale(int nrow, int ncol, float **ma, float sFactor);
void matrix_shift(int nrow, int ncol, float **ma, float sFactor);
void matrix_sum(int nrow, int ncol, float **ma, float **ma2);
void parse_comm_args(int *argc, char *argv[], t_mdPara *md_para,           
                    char *fn_in, char *fn_out);
void write_xyz(FILE *fp, t_mdPara *md_para, float **x);
void write_all(FILE *fp, t_mdPara *md_para, t_state *md_state);
void read_xyz(char *fn, t_state *md_state, t_mdPara *md_para);
int read_traj(FILE *fp, t_state *md_state, t_mdPara *md_para);
void get_vec_pbc(float length, fvec x1, fvec x2, fvec vecij);
void gen_velocity(t_mdPara *md_para, t_state *md_state);
void calc_kinE(t_mdPara *md_para, t_state *md_state);
float calc_force(t_mdPara *md_para, fvec *x, fvec *f);
void update(t_mdPara *md_para, t_state *md_state);
