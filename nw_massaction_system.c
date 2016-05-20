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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_matrix.h>

#include "nw_massaction.h"
#include "nw_massaction_system.h"

#define OUT stdout
#define ERR stderr

/* Allocate memory for ncomp compounds and nreac reactions */
nw_massaction_system *nw_massaction_system_alloc(int ncomp, int nreact, 
						 const char* names[]) {
  nw_massaction_system* sys=malloc(sizeof(*sys));

  sys->ncomp=ncomp;
  sys->nreact=nreact;

  sys->last=0;

  sys->names=malloc(sys->ncomp*sizeof(*(sys->names)));
  sys->react=malloc(sys->nreact*sizeof(*(sys->react)));

  nw_massaction_system_setNames(sys, names,ncomp);

  return sys;
}

/* Returns number of compounds */
size_t nw_massaction_system_ncomp(const nw_massaction_system *sys) {
  return sys->ncomp;
}

/* Returns number of reaction */
size_t nw_massaction_system_nreact(const nw_massaction_system *sys) {
  return sys->last;
}

/* Set k-th reaction in mass action system */
void nw_massaction_system_setReact(nw_massaction_system *sys, int k, nw_massaction *r) {
  sys->react[k]=r;
}

/* Set k-th reaction in mass action system */
nw_massaction *nw_massaction_system_getReact(const nw_massaction_system *sys, int k) {
  return sys->react[k];
}

/* Set name of k-th compound */
void nw_massaction_system_setName(nw_massaction_system *sys, int k, const char *name) {
  sys->names[k]=name;
}

/* Set name of k-th compound */
void nw_massaction_system_setNames(nw_massaction_system *sys, const char *name[], size_t N) {
  int i;
  for(i=0; i<N; i++) nw_massaction_system_setName(sys, i, name[i]);
}

/* Get name of k-th compound */
char * nw_massaction_system_getName(const nw_massaction_system *sys, int k) {
  return sys->names[k];
}

/* Adds reaction to nw_massaction_system at index k0 */
void nw_massaction_system_addReactK0(nw_massaction_system *sys, 
				   int k0,
				   int i, int j, int k, 
				   double kF, double kB) {
  nw_massaction *r;
  char *nameA=nw_massaction_system_getName(sys, i);
  char *nameB=nw_massaction_system_getName(sys, j);
  char *nameC=nw_massaction_system_getName(sys, k);

  r=nw_massaction_ABC(i, j, k, nameA, nameB, nameC, kF, kB);

  nw_massaction_system_setReact(sys, k0, r);
}

/* Adds reaction to nw_massaction_system at index k0 */
void nw_massaction_system_addEnzymeK0(nw_massaction_system *sys, 
				      int k0,
				      int i, int j, int k, 
				      double kF) {
  nw_massaction *r;
  char *nameA=nw_massaction_system_getName(sys, i);
  char *nameB=nw_massaction_system_getName(sys, j);
  char *nameC=nw_massaction_system_getName(sys, k);


  r=nw_massaction_enzyme(i, j, k, nameA, nameB, nameC, kF);
  nw_massaction_system_setReact(sys, k0, r);
}

/* Adds reaction to nw_massaction_system */
void nw_massaction_system_addReact(nw_massaction_system *sys, 
				   int i, int j, int k, 
				   double kF, double kB) {
  double k0=sys->last;

  nw_massaction_system_addReactK0(sys, k0, i, j, k, kF, kB);
  sys->last++;
}

/* Adds enzymatic reaction to nw_massaction_system */
void nw_massaction_system_addEnzyme(nw_massaction_system *sys, 
				    int i, int j, int k, 
				    double kF) {
  double k0=sys->last;

  nw_massaction_system_addEnzymeK0(sys, k0, i, j, k, kF);
  sys->last++;
}

/* Adds reaction to nw_massaction_system where A is the complex on the
   left-hand-side and B is the complex of the right-hand side of the
   equation. */
void nw_massaction_system_addGeneralReaction(nw_massaction_system *sys, 
					     const int* Ak, const int *AS,
					     size_t nA,
					     const int *Bk, const int *BS,
					     size_t nB,
					     double kF, double kB) {
  int i;
  double k0=sys->last;
  nw_massaction *ma=nw_massaction_alloc(nA, nB, kF, kB);

  for(i=0; i<nA; i++) {
    char *name=nw_massaction_system_getName(sys, Ak[i]);
    /* Set reactant i in complex A */
    nw_massaction_setA(ma, i, name, Ak[i], AS[i]);
  }

  for(i=0; i<nB; i++) {
    char *name=nw_massaction_system_getName(sys, Bk[i]);
    /* Set reactant i in complex B */
    nw_massaction_setB(ma, i, name, Bk[i], BS[i]);
  }

  nw_massaction_system_setReact(sys, k0, ma);
  sys->last++;
}

/* Adds reaction to nw_massaction_system where A is the complex on the
   left-hand-side and B is the complex of the right-hand side of the
   equation. */
void nw_massaction_system_addGeneralReaction_S1(nw_massaction_system *sys, 
						const int* Ak, size_t nA,
						const int *Bk, size_t nB,
						double kF, double kB) {
  int i;
  double k0=sys->last;
  nw_massaction *ma=nw_massaction_alloc(nA, nB, kF, kB);

  for(i=0; i<nA; i++) {
    char *name=nw_massaction_system_getName(sys, Ak[i]);
    /* Set reactant i in complex A */
    nw_massaction_setA(ma, i, name, Ak[i], 1);
  }

  for(i=0; i<nB; i++) {
    char *name=nw_massaction_system_getName(sys, Bk[i]);
    /* Set reactant i in complex B */
    nw_massaction_setB(ma, i, name, Bk[i], 1);
  }

  nw_massaction_system_setReact(sys, k0, ma);
  sys->last++;
}

/* Updates kF in reaction i of nw_massaction_system */
void nw_massaction_system_setReact_kF(nw_massaction_system *sys, 
				      int i, double kF) {
  nw_massaction *r=nw_massaction_system_getReact(sys, i);

  nw_massaction_set_kF(r, kF);
}

/* Updates kB in reaction i of nw_massaction_system */
void nw_massaction_system_setReact_kB(nw_massaction_system *sys, 
				      int i, double kB) {
 nw_massaction *r=nw_massaction_system_getReact(sys, i);

  nw_massaction_set_kB(r, kB);
}

/* Print chemical reactions of mass action system */
void nw_massaction_system_print(FILE *fp, const nw_massaction_system * sys) {
  int k;
  for(k=0; k<sys->last; k++) {
    nw_massaction *r=nw_massaction_system_getReact(sys, k);

    nw_massaction_print(fp, r);
  }
}

/* Returns stoichiometric matrix */
gsl_matrix *nw_massaction_system_S(const nw_massaction_system * sys) {
  int i, j;
  int nR=nw_massaction_system_nreact(sys);
  gsl_matrix *S=gsl_matrix_calloc(nw_massaction_system_ncomp(sys), 
				  nR);

  for(i=0; i<nR; i++) {
    nw_massaction *r=nw_massaction_system_getReact(sys, i);
    /* Complex A*/
    for(j=0; j<r->nA; j++) {
      int rInd=nw_massaction_getAr(r, j);
      int s=nw_massaction_getAs(r, j);
      gsl_matrix_set(S, rInd, i, s);
    }

    /* Complex B*/
    for(j=0; j<r->nB; j++) {
      int rInd=nw_massaction_getBr(r, j);
      int s=nw_massaction_getBs(r, j);
      int Sri=gsl_matrix_get(S, rInd, i);
      gsl_matrix_set(S, rInd, i, Sri-s);
    }
  }
  
  return S;
}

/* Returns derivative dydt for state vector y where par is a nw_massaction_system */
int nw_massaction_system_ODEs(double t, const double y[], double dydt[], void *par) {
  int i;
  nw_massaction_system *sys=(nw_massaction_system *)par;
  size_t nR=sys->last;
  size_t nC=sys->ncomp;

  for(i=0; i<nC; i++) {
    dydt[i]=0;
  }

  for(i=0; i<nR; i++) {
    nw_massaction *r=nw_massaction_system_getReact(sys, i);
    nw_massaction_setODE(r, y, dydt);
  }
  
  return GSL_SUCCESS;
}

/* Allocate space for solving full system */
void nw_massaction_system_setupReducedODEs(nw_massaction_system *sys) {
  size_t nC=nw_massaction_system_ncomp(sys);
  sys->yFull=malloc(nC*sizeof(*(sys->yFull)));
  sys->dydtFull=malloc(nC*sizeof(*(sys->dydtFull)));
}

/* Allocate space for solving full system */
void nw_massaction_system_setupJac(nw_massaction_system *sys) {
  size_t nC=nw_massaction_system_ncomp(sys);
  sys->J=malloc(nC*nC*sizeof(*(sys->J)));
}

/* Add conserved quantities. N should be nComp */
void nw_massaction_system_addCons(nw_massaction_system *sys, const double *kernelS, int M, int N) {
  sys->cons = nw_cons_init(kernelS, M, N);
  nw_massaction_system_setupReducedODEs(sys);
  nw_massaction_system_setupJac(sys);
}

nw_cons *nw_massaction_system_getConserved(const nw_massaction_system *sys) {
  return sys->cons;
}

/* Initialise conserved quantities from state vector y */
double * nw_massaction_system_conservedQuantities(nw_massaction_system *sys, const double *y) {
  nw_cons *cons =nw_massaction_system_getConserved(sys);
  if (cons)
    return nw_cons_conservedQuantities(cons, y);
  return NULL;
}

/* Set conserved quantity i to new value cInew */
void nw_massaction_system_setCons(nw_massaction_system *sys, int i, double cInew) {
  nw_cons *c=nw_massaction_system_getConserved(sys);
  nw_cons_setCons(c, i, cInew);
}

/* Get conserved quantity i */
double nw_massaction_system_getCons(const nw_massaction_system *sys, int i) {
  nw_cons *c=nw_massaction_system_getConserved(sys);
  return nw_cons_getCons(c, i);
}

/* Calculate reduced system of par with conserved quantities removed. It is assumed that y and dydt have allocated sufficient */
int nw_massaction_system_reducedODEs(double t,
			const double y[],
			double dydt[], void *par) {
  nw_massaction_system *sys=(nw_massaction_system *)par;
  int err;
  nw_cons *consABC=nw_massaction_system_getConserved(sys);

  /* Initialise vector of full system with independent quantities. */
  /* for(i=0; i<nC-nCons; i++) { */
  /*   yFull[i+nCons]=y[i]; */
  /* } */
  nw_cons_copyIndependentToFull(consABC, sys->yFull, y);
  /* copyIndependentToFull(&Gview.matrix, yFull, y); */
  /* fprintf(OUT, "nw_cons_reducedODES: Copied independent quantities...\n"); */

  nw_cons_xCons(consABC, sys->yFull);
  /* fprintf(OUT, "nw_cons_reducedODES: Calculated dependent quantities...\n"); */

  err=nw_massaction_system_ODEs(t, sys->yFull, sys->dydtFull, par);
  /* fprintf(OUT, "nw_cons_reducedODES: Solved full system...\n"); */


  if(err != GSL_SUCCESS) return err;

  /* for(i=0; i<nC-nCons; i++) { */
  /*   dydt[i]=dydtFull[i+nCons]; */
  /* } */
  nw_cons_copyIndependentFromFull(consABC, dydt, sys->dydtFull);
  /* fprintf(OUT, "nw_cons_reducedODES: Copied independent derivatives...\n"); */
  return GSL_SUCCESS;
}

/* Copy 'independent quantities' (according to conservation laws) from
   vector vFull to v (essentially copies indices that are not lead
   coefficients). */
void nw_massaction_system_copyIndependentFromFull(const nw_massaction_system *sys, 
						  double *v, const double *vFull) {
  nw_cons *c=nw_massaction_system_getConserved(sys);
  nw_cons_copyIndependentFromFull(c, v, vFull);
}

/* Copy 'independent quantities' (according to conservation laws) from
   vector vFull to v (essentially copies indices that are not lead
   coefficients). */
void nw_massaction_system_copyIndependentToFull(const nw_massaction_system *sys, 
						double *vFull, const double *v) {
  nw_cons *c=nw_massaction_system_getConserved(sys);
  nw_cons_copyIndependentToFull(c, vFull, v);
}

/* Print names of reactants in massaction_system sys */
void nw_massaction_system_printNames(FILE *fp, const nw_massaction_system *sys) {
  int i, nC=nw_massaction_system_ncomp(sys);
  for(i=0; i<nC-1; i++) 
    fprintf(fp, "%s\t", nw_massaction_system_getName(sys, i));
  fprintf(fp, "%s\n", nw_massaction_system_getName(sys, i));
}

/* Print names of 'independent quantities' */
void nw_massaction_system_printIndependentNames(FILE *fp, const nw_massaction_system *sys) {
  nw_cons *c=nw_massaction_system_getConserved(sys);
  int rows=nw_cons_getNcons(c), cols=nw_cons_getNreact(c);
  int i, j;
  int lead=nw_cons_getLead(c, 0), nextLead;
  for(i=1; i<=rows; i++) {
    /*Print everything between lead coefficients, last "lead coefficient" is the last column*/
    if(i==rows) nextLead=cols;
    else nextLead=nw_cons_getLead(c, i);
    for(j=lead+1; j<nextLead; j++) {
      /* fprintf(fp, "%s(%i)\t", nw_massaction_system_getName(sys, j), j); */
      fprintf(fp, "%s\t", nw_massaction_system_getName(sys, j));
    } 
    lead=nextLead;
  }
  fprintf(fp, "\n");
}

/* Print names of 'dependent quantities' */
void nw_massaction_system_printDependentNames(FILE *fp, const nw_massaction_system *sys) {
  nw_cons *c=nw_massaction_system_getConserved(sys);
  int rows=nw_cons_getNcons(c);
  int i;
  for(i=0; i<rows; i++) {
    int lead=nw_cons_getLead(c, i);
    /* fprintf(fp, "%s(%i)\t", nw_massaction_system_getName(sys, lead), lead); */
    fprintf(fp, "%s\t", nw_massaction_system_getName(sys, lead));
  } 
  fprintf(fp, "\n");
}

/* Set Jacobian for reactions of massaction_system sys. */
void nw_massaction_system_setJac(const nw_massaction_system *sys, const double *y, double *dfdy) {
  int i;
  size_t nC=nw_massaction_system_ncomp(sys), nR=nw_massaction_system_nreact(sys);

  for(i=0; i<nC*nC; i++) {
    dfdy[i]=0;
  }
  for(i=0; i<nR; i++) {
    nw_massaction *r=nw_massaction_system_getReact(sys, i);
    nw_massaction_setJac(r, y, dfdy, nC);
  }
}

/* Save Jacobian at point y in dfdy. Parameter par contains
   nw_massaction_sys.  */
int nw_massaction_system_jac(double t, const double *y, double *dfdy, 
			     double dfdt[], void *par) {
  const nw_massaction_system *sys=(const nw_massaction_system *)par;
  int i;

  nw_massaction_system_setJac(sys, y, dfdy);
  for(i=0; i<nw_massaction_system_nreact(sys); i++) dfdt[i]=0;
  return GSL_SUCCESS;
}

/* Calculates total derivative of f_i by x_j for with all dependent quantities removed. */
double nw_massaction_system_totalDerivative(const nw_massaction_system *sys, 
				     const double *dfdy, 
				     int i, int j) {
  nw_cons *cons=nw_massaction_system_getConserved(sys);
  int nC=nw_massaction_system_ncomp(sys);
  int nCons=nw_cons_getNcons(cons);
  gsl_matrix_const_view dfdy_mat= gsl_matrix_const_view_array (dfdy, nC, nC);
  const gsl_matrix * J = &dfdy_mat.matrix;

  double Jij=gsl_matrix_get(J, i, j);
  int k;
  /* Derivatives by dependent quantities */
  for(k=0; k<nCons; k++) {
    int depK=nw_cons_getLead(cons, k);
    /* fprintf(ERR, "nw_massaction_system_totalDerivative(): k=%i, depK=%i\n", k, depK); */
    double dXdk_dXj=nw_cons_dX_Li_d_Xj(cons, k, j);
    Jij+=dXdk_dXj*gsl_matrix_get(J, i, depK);
    /* fprintf(ERR,"nw_massaction_system_totalDerivative(): J(%i,%i)=%g\n", depK, j, Jij); */
  }
  return Jij;
}

/* Set Jacobian for reduced reactions of massaction_system sys. */
void nw_massaction_system_setReducedJac(const nw_massaction_system *sys, const double *y, double *dfdy) {
  nw_cons *cons=nw_massaction_system_getConserved(sys);
  int i,j;
  int nC=nw_massaction_system_ncomp(sys);
  int nCons=nw_cons_getNcons(cons);
  /* gsl_matrix_view J_mat= gsl_matrix_view_array (sys->J, nC, nC); */
  /* gsl_matrix * J = &J_mat.matrix; */
  gsl_matrix_view dfdy_mat= gsl_matrix_view_array (dfdy, nC-nCons, nC-nCons);
  gsl_matrix * Jred = &dfdy_mat.matrix;

  for(i=0; i<nC; i++) sys->yFull[i]=0;
  for(i=0; i<(nC-nCons)*(nC-nCons); i++) dfdy[i]=0;

  /* Calculate FULL vector of concentrations from conserved quantities. */
  nw_cons_copyIndependentToFull(cons, sys->yFull, y);
  /* copyIndependentToFull(&Gview.matrix, yFull, y); */
  /* fprintf(OUT, "nw_cons_reducedJac: Copied independent quantities...\n"); */
  nw_cons_xCons(cons, sys->yFull);
  /* for(i=0; i<nC; i++) { */
  /*   fprintf(ERR, "%g\t", sys->yFull[i]); */
  /* } */
  /* fprintf(ERR, "\n"); */

  /* Get Jacobian for full system */
  /* fprintf(OUT, "nw_cons_reducedJac: Calculate full Jacobian...\n"); */
  nw_massaction_system_setJac(sys, sys->yFull, sys->J);

  /* fprintf(OUT, "nw_cons_reducedJac(): Calculate reduced Jacobian...\n"); */
  /* fprintf(ERR, "Lead coefficients: "); nw_cons_printLeads(ERR, cons); */
  /* fprintf(ERR, "Free coefficients: "); nw_cons_printFree(ERR, cons); */

  for(i=0; i<nC-nCons; i++) {
      for(j=0; j<nC-nCons; j++) {
	int indI=nw_cons_getFree(cons, i);
	int indJ=nw_cons_getFree(cons, j);
	double Jij=nw_massaction_system_totalDerivative(sys, sys->J, indI, indJ);
	/* fprintf(ERR, "nw_cons_reducedJac(): (%i,%i): Jred(%i,%i)=%g\n", i, j, indI, indJ, Jij); */
	gsl_matrix_set(Jred, i, j, Jij);
      }    
  }
}

/* Save reduced Jacobian at point y in dfdy. Parameter par contains
   nw_massaction_sys.  */
int nw_massaction_system_reducedJac(double t, const double *y, double *dfdy, 
				    double dfdt[], void *par) {
  const nw_massaction_system *sys=(const nw_massaction_system *)par;
  int i;
  nw_massaction_system_setReducedJac(sys, y, dfdy);
  for(i=0; i<nw_massaction_system_nreact(sys); i++) dfdt[i]=0;
  return GSL_SUCCESS;
}
