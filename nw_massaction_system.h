#ifndef NW_MA_SYS
#define NW_MA_SYS
/*
 * nw_massaction_system.c
 *
 *
 * Simple mass action kinetics
 *
 * Ivo Siekmann, 07/03/2014
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <gsl/gsl_matrix.h>
#include "nw_massaction.h"
#include "nw_cons.h"

typedef struct {
  /* Number of compounds */
  size_t ncomp;
  /* Number of reactions */
  size_t nreact;

  /* Index of last reaction */
  size_t last;
  
  /* Names of compounds (ncomp)*/
  char **names;
  /* Reactions (nReac) */
  nw_massaction ** react;

  /* Conserved quantities */
  nw_cons *cons;

  /* State vector and derivatives of full equations */
  double *yFull, *dydtFull;

  /*Full Jacobian */
  double *J;
} nw_massaction_system;

/* Allocate memory for ncomp compounds and nreac reactions */
nw_massaction_system *nw_massaction_system_alloc(int ncomp, int nreact, const char *names[]);

/* Returns number of compounds */
size_t nw_massaction_system_ncomp(const nw_massaction_system *sys);

/* Returns number of reaction */
size_t nw_massaction_system_nreact(const nw_massaction_system *sys);

/* Set k-th reaction in mass action system */
void nw_massaction_system_setReact(nw_massaction_system *sys, int k, nw_massaction *r);

/* Returns derivative dydt for state vector y where par is a nw_massaction_system */
int nw_massaction_system_ODEs(double t, const double y[], double dydt[], void *par);

/* Set k-th reaction in mass action system */
nw_massaction *nw_massaction_system_getReact(const nw_massaction_system *sys, int k);

/* Set name of k-th compound */
void nw_massaction_system_setName(nw_massaction_system *sys, int k, const char *name);

/* Set name of k-th compound */
void nw_massaction_system_setNames(nw_massaction_system *sys, const char *name[], size_t N);

/* Get name of k-th compound */
char * nw_massaction_system_getName(const nw_massaction_system *sys, int k);

/* Adds reaction to nw_massaction_system */
void nw_massaction_system_addReact(nw_massaction_system *sys, 
				   int i, int j, int k, 
				   double kF, double kB);

/* Adds reaction to nw_massaction_system where A is the complex on the
   left-hand-side and B is the complex of the right-hand side of the
   equation. */
void nw_massaction_system_addGeneralReaction(nw_massaction_system *sys, 
					     const int* Ak, const int *AS,
					     size_t nA,
					     const int *Bk, const int *BS,
					     size_t nB,
					     double kF, double kB);

/* Adds reaction to nw_massaction_system where A is the complex on the
   left-hand-side and B is the complex of the right-hand side of the
   equation. */
void nw_massaction_system_addGeneralReaction_S1(nw_massaction_system *sys, 
						const int* Ak, size_t nA,
						const int *Bk, size_t nB,
						double kF, double kB);

/* Updates kF in reaction i of nw_massaction_system */
void nw_massaction_system_setReact_kF(nw_massaction_system *sys, 
				      int i, double kF);

/* Updates kB in reaction i of nw_massaction_system */
void nw_massaction_system_setReact_kB(nw_massaction_system *sys, 
				      int i, double kB);

/* Adds enzymatic reaction to nw_massaction_system */
void nw_massaction_system_addEnzyme(nw_massaction_system *sys, 
				    int i, int j, int k, 
				    double kF);

/* Print chemical reactions of mass action system */
void nw_massaction_system_print(FILE *fp, const nw_massaction_system * sys);

/* Returns stoichiometric matrix */
gsl_matrix *nw_massaction_system_S(const nw_massaction_system * sys);

/* Returns derivative dydt for state vector y where par is a nw_massaction_system */
int nw_massaction_system_ODEs(double t, const double y[], double dydt[], void *par);

/* Add conserved quantities */
void nw_massaction_system_addCons(nw_massaction_system *sys, const double *kernelS, int M, int N);

/* Set conserved quantity i to new value cInew */
void nw_massaction_system_setCons(nw_massaction_system *sys, int i, double cInew);

/* Get conserved quantity i */
double nw_massaction_system_getCons(const nw_massaction_system *sys, int i);

/* Initialise conserved quantities from state vector y */
/* void nw_massaction_system_conservedQuantities(nw_massaction_system *sys, const double *y); */
double * nw_massaction_system_conservedQuantities(nw_massaction_system *sys, const double *y);

/* Calculate reduced system of par with conserved quantities removed. It is assumed that y and dydt have allocated sufficient */
int nw_massaction_system_reducedODEs(double t, const double y[], double dydt[], void *par);

/* Copy 'independent quantities' (according to conservation laws) from
   vector vFull to v (essentially copies indices that are not lead
   coefficients). */
void nw_massaction_system_copyIndependentFromFull(const nw_massaction_system *sys, 
						  double *v, const double *vFull);

/* Copy 'independent quantities' (according to conservation laws) from
   vector vFull to v (essentially copies indices that are not lead
   coefficients). */
void nw_massaction_system_copyIndependentToFull(const nw_massaction_system *sys, 
						double *vFull, const double *v);

/* Print names of reactants in massaction_system sys */
void nw_massaction_system_printNames(FILE *fp, const nw_massaction_system *sys);

/* Print names of 'independent quantities' */
void nw_massaction_system_printIndependentNames(FILE *fp, const nw_massaction_system *sys);

/* Print names of 'dependent quantities' */
void nw_massaction_system_printDependentNames(FILE *fp, const nw_massaction_system *sys);

/* Save Jacobian at point y in dfdy. Parameter par contains
   nw_massaction_sys.  */
int nw_massaction_system_jac(double t, const double *y, double *dfdy, double dfdt[],  void *par);

/* Calculates total derivative of f_i by x_j for with all dependent quantities removed. */
double nw_massaction_system_totalDerivative(const nw_massaction_system *sys, const double *dfdy, 
					    int i, int j);

/* Set Jacobian for reduced reactions of massaction_system sys. */
void nw_massaction_system_setReducedJac(const nw_massaction_system *sys, const double *y, double *dfdy);

/* Save reduced Jacobian at point y in dfdy. Parameter par contains
   nw_massaction_sys.  */
int nw_massaction_system_reducedJac(double t, const double *y, double *dfdy, double dfdt[], void *par);
#endif
