#ifndef NW_CONS
#define NW_CONS
/*
 * nw_cons.h
 *
 * Conserved quantities in a mass action kinetics based model.
 * 
 * 
 *
 * Ivo Siekmann, 23.05.2014
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

typedef struct {
  /* Indices of lead coefficients for each 'row' */
  int *leads;

  /* Indices of non-lead coefficients = independent (free) quantities */
  int *free;

  /* kernel of the stoichiometric matrix S in row echelon form contains
   all information on conserved quantities */
  double **kerS;
  size_t rows,cols;

  /* Vector of conserved quantities */
  double *Cons;
} nw_cons;

/* Returns the number of conserved quantities */
int nw_cons_getNcons(const nw_cons *c);

/* Returns the number of reactants */
int nw_cons_getNreact(const nw_cons *c);

void nw_cons_setKernelS(nw_cons *c, const double *kernelS, int M, int N);

void nw_cons_initCons(nw_cons *c);

/* Allocates nw_cons instance Conserved quantities are initialised by nw_cons_conservedQuantities*/
nw_cons* nw_cons_init(const double *kernelS, int M, int N);

/* Set conserved quantity i to C */
void nw_cons_setCons(nw_cons *c, int i, double C);

/* Get conserved quantity i */
double nw_cons_getCons(const nw_cons *c, int i);

/*Get column index of lead coefficient in row i of matrix R. Assuming
  R is in row echelon form. */
int nw_cons_getLead(const nw_cons *c, int i);

/*Get column index of lead coefficient in row i of matrix R. */
int nw_cons_getFree(const nw_cons *c, int i);

/* Calculate conserved quantities from state vector y and initialise
   in nw_cons. Return the result. */
double *nw_cons_conservedQuantities(nw_cons *c, const double *y);

/* Calculate concentration of n-th DEPENDENT reactant from conserved
   quantity T. n must be smaller than number of conservation laws. */
double nw_cons_xNfromCons(const nw_cons *c, const double* x, int n);

/*Calculate all conserved quantities by BACK SUBSTITUTION. The
  INDEPENDENT quantities set in the vector x are substituted, the
  dependent quantities are ignored and replaced by what is calculated
  from the conservation laws. */
void nw_cons_xCons(const nw_cons *c, double* x);

/* Returns derivative of DEPENDENT quantity i by independent quantity j */
double nw_cons_dX_Li_d_Xj(const nw_cons *c, int i, int j);

/* Copy 'independent quantities' (according to conservation laws) from
   vector vFull to v (essentially copies indices that are not lead
   coefficients). */
void nw_cons_copyIndependentFromFull(const nw_cons *c, double *v, const double *vFull);

/* Copy 'independent quantities' (according to conservation laws) from
   vector vFull to v (essentially copies indices that are not lead
   coefficients). */
void nw_cons_copyIndependentToFull(const nw_cons *c, double *vFull, const double *v);

/* Deallocate memory */
void nw_cons_free(nw_cons *c);

/* Calculate reduced system of par with conserved quantities c
   removed. It is assumed that y and dydt have allocated sufficient */
/* int nw_cons_reducedODEs(const nw_cons *c, */
/* 			double t, */
/* 			const double y[], */
/* 			double dydt[], void *par); */
#endif
 
