/*
 * nw_cons.c
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

#include <stdio.h>
#include <stdlib.h>

#define OUT stdout
#define ERR stderr


#include <gsl/gsl_errno.h>

#include "nw_massaction_system.h"
#include "nw_cons.h"

/* kernel of the stoichiometric matrix S in row echelon form contains
   all information on conserved quantities */
/* static double **kerS;/\* ={1,0,1, *\/ */
/* /\* 	0,1,1}; *\/ */
/* static int rows, cols; */

/* /\* Two conserved quantities *\/ */
/* static double *Cons; */

/* For solving ODE systems */
static double *yFull, *dydtFull;

/* Returns the number of conserved quantities */
int nw_cons_getNcons(const nw_cons *c) {
  return c->rows;
}

/* Returns the number of reactants */
int nw_cons_getNreact(const nw_cons *c) {
  return c->cols;
}

/* Get kernel component Kij */
double nw_cons_getKij(const nw_cons *c, int i, int j) {
  return c->kerS[i][j];
}


/*Get column index of lead coefficient in row i of matrix R. */
int nw_cons_getLead(const nw_cons *c, int i) {
  return c->leads[i];
}

/*Set lead coefficient i to l.*/
void nw_cons_setLead(nw_cons *c, int i, int l) {
  c->leads[i]=l;
}

/*Get column index of lead coefficient in row i of matrix R. */
int nw_cons_getFree(const nw_cons *c, int i) {
  return c->free[i];
}

/*Set lead coefficient i to l.*/
void nw_cons_setFree(nw_cons *c, int i, int f) {
  c->free[i]=f;
}

/* Calculate index of lead coefficient in row i of matrix R. Assuming
  R is in row echelon form. */
int nw_cons_calculateLead(const nw_cons *c, int i) {
  int j=i;
  double rIJ=nw_cons_getKij(c, i, j);
  if(!(rIJ==0)) return j;
  /* Find first non-zero entry in row i starting from the diagonal */
  for(j=i; /* fprintf(OUT, "nw_cons_getLead(): R(%i,%i)=%f\n", i, j, rIJ), */ (rIJ == 0) && (j<nw_cons_getNreact(c)); j++) {
    rIJ=nw_cons_getKij(c, i, j);
  }
  return j-1;
}

/* Initialise vector of lead coefficients */
void nw_cons_initLeads(nw_cons *c) {
  int rows=nw_cons_getNcons(c);
  int i;
  
  for(i=0; i<rows; i++) {
    int leadI=nw_cons_calculateLead(c, i);
    nw_cons_setLead(c, i, leadI);
  }
}

/* Initialise vector of "free" coefficients (assumes that lead coefficients are already calculated!) */
void nw_cons_initFree(nw_cons *c) {
  int rows=nw_cons_getNcons(c), cols=nw_cons_getNreact(c);
  int i, j, k=0;
  int lead=nw_cons_getLead(c, 0), nextLead;
  
  for(i=1; i<=rows; i++) {
    /*Print everything between lead coefficients, last "lead coefficient" is the last column*/
    if(i==rows) nextLead=cols;
    else nextLead=nw_cons_getLead(c, i);
    for(j=lead+1; j<nextLead; j++) {
      /* fprintf(stderr, "lead=%i, nextLead=%i, Free(%i)=%i\t", lead, nextLead, k, j); */
      nw_cons_setFree(c, k++, j);
    } 
    lead=nextLead;
  }
  /* fprintf(stderr, "\n"); */
}

/* Print indices of lead coefficients */
void nw_cons_printLeads(FILE *fp, const nw_cons *c) {
  int i;
  int nC=nw_cons_getNcons(c);
  int lI;
  for(i=0; i<nC-1; i++) {
    lI=nw_cons_getLead(c, i);
    fprintf(fp, "%i\t", lI);
  }
  lI=nw_cons_getLead(c, i);
  fprintf(fp, "%i\n", lI);
}

/* Print indices of free coefficients */
void nw_cons_printFree(FILE *fp, const nw_cons *c) {
  int i;
  int nF=nw_cons_getNreact(c)-nw_cons_getNcons(c);
  int fI;
  for(i=0; i<nF-1; i++) {
    fI=nw_cons_getFree(c, i);
    fprintf(fp, "%i\t", fI);
  }
  fI=nw_cons_getFree(c, i);
  fprintf(fp, "%i\n", fI);
}


void nw_cons_setKernelS(nw_cons *c, const double *kernelS, int M, int N) {
  int i, j;
  for(i=0; i<c->rows; i++) {
    for(j=0; j<c->cols; j++) {
      c->kerS[i][j]=kernelS[i*N+j];
    }
  }
}

void nw_cons_initCons(nw_cons *c) {
  size_t N=nw_cons_getNcons(c);

  c->Cons=malloc(N*sizeof(*(c->Cons)));
  /* c->Cons=malloc(N*sizeof(double)); */
}

/* Just leave in temporarily! Should be done in massaction system! */
void nw_cons_allocODE(nw_cons *c) {
  size_t N=nw_cons_getNreact(c);

  yFull=malloc(N*sizeof(*yFull));
  dydtFull=malloc(N*sizeof(*dydtFull));
}


nw_cons *nw_cons_init(const double *kernelS, int M, int N) {
  nw_cons *c=malloc(sizeof(*c));
  int i;
  c->rows=M, c->cols=N;

  /* Allocate memory for lead coefficients */
  c->leads=malloc(c->rows*sizeof(*(c->leads)));
  /* Allocate memory for free coefficients */
  c->free=malloc((N-M)*sizeof(*(c->free)));


  /* Allocate memory for kernel matrix */
  c->kerS=malloc(c->rows*sizeof(*(c->kerS)));
  for(i=0; i<c->rows; i++) {
    c->kerS[i]=malloc(c->cols*sizeof(**(c->kerS)));
  }

  /* Allocate memory for conserved quantities */
  nw_cons_initCons(c);

  /* Allocate memory for ODE solving without reduced quantities */
  /* nw_cons_allocODE(); */

  
  nw_cons_setKernelS(c, kernelS, M, N);

  /*Initialise lead coefficients - must be done after kernel is set!*/
  /* fprintf(stderr, "nw_cons_init(): Initialising dependent quantities:\n"); */
  nw_cons_initLeads(c);
  /* fprintf(stderr, "nw_cons_init(): Initialised dependent quantities:\n"); */
  /* nw_cons_printLeads(stderr, c); */


  /*Initialise free coefficients - must be done lead coefficients are set!*/
  /* fprintf(stderr, "nw_cons_init(): Initialising independent quantities:\n"); */
  nw_cons_initFree(c);

  /* fprintf(stderr, "nw_cons_init(): Initialised independent quantities:\n"); */
  /* nw_cons_printFree(stderr, c); */

  return c;
}

/* Deallocate memory */
void nw_cons_free(nw_cons *c) {
  int i;
  size_t rows=nw_cons_getNcons(c);

  free(c->leads);
  free(c->free);

  for(i=0; i<rows; i++) {
    free(c->kerS[i]);
  }
  free(c->kerS);

  free(c->Cons);
}

/* Set conserved quantity i to C */
void nw_cons_setCons(nw_cons *c, int i, double C) {
  c->Cons[i]=C;
}

/* Get conserved quantity i */
double nw_cons_getCons(const nw_cons *c, int i) {
  double out=c->Cons[i];
  /* fprintf(OUT, "nw_cons_getCons(): T(%i)=%g\n", i, out); */
  return out;
}


/* Calculate conserved quantities from state vector y. Returns the result. */
double *nw_cons_conservedQuantities(nw_cons *c, const double *y) {
  int i, j;

  for(i=0; i<c->rows; i++) {
    double cI=0;
    for(j=0; j<c->cols; j++) {
      cI+=nw_cons_getKij(c, i,j)*y[j];
    }
    fprintf(OUT, "nw_cons_conservedQuantities(): T(%i) = %g\n", i, cI);
    nw_cons_setCons(c, i, cI);
  }

  return c->Cons;
}

/* Calculate concentration of n-th DEPENDENT reactant from conserved
   quantity T. n must be smaller than number of conservation laws. */
double nw_cons_xNfromCons(const nw_cons *c, const double* x, int n) {
  /* size_t k=G->size1, m=G->size2; */
  size_t k=nw_cons_getLead(c, n)+1, m=nw_cons_getNreact(c);
  int j;
  double xC;/* =nw_cons_getCons(c, n); */

  /* fprintf(OUT, "nw_cons_xNfromCons(): T(%i)=", n); */
  xC=nw_cons_getCons(c, n);
  /* fprintf(OUT,"%g\n", xC); */

  /* Subtract everything right of index k. It is assumed that coefficient of index k is 1! */
  for(j=k; j<m; j++) {
    /* fprintf(OUT, "nw_cons_xNfromCons(): ker(%i,%i)=%g, x[%i]=%g\n", n, j, nw_cons_getKij(c, n, j), j, x[j]); */
    xC -= nw_cons_getKij(c, n, j)*x[j];
  }
  return xC;
}

/*Calculate all conserved quantities by BACK SUBSTITUTION. The
  INDEPENDENT quantities set in the vector x are substituted, the
  dependent quantities are ignored and replaced by what is calculated
  from the conservation laws. */
void nw_cons_xCons(const nw_cons *c, double* x) {
  int k=nw_cons_getNcons(c);
  int i;
  /* fprintf(OUT, "xCons():\n"); */
  for(i=0; i<k; i++) {
    int m=nw_cons_getLead(c, k-1-i);
    double xI;
    /* fprintf(OUT, "xCons(): %i-th lead coefficient: %i\n", k-1-i, m); */
    xI=nw_cons_xNfromCons(c, x, k-1-i);
    /* Dependent quantity is where lead coefficient is */
    x[m]=xI;
  }
}

/* Returns index of next independent quantity after i */
/* int nw_cons_nextIndependent(const nw_cons *c, int i) { */
/*   int rows=nw_cons_getNcons(c), cols=nw_cons_getNreact(c); */
/*   int i, j, k=0; */
/*   int lead=nw_cons_getLead(c, 0), nextLead; */
/*   for(i=1; i<=rows; i++) { */
/*     /\*Copy everything between lead coefficient, last "lead coefficient" is the last column*\/ */
/*     if(i==rows) nextLead=cols; */
/*     else nextLead=nw_cons_getLead(c, i); */
/* } */

/* Copy 'independent quantities' (according to conservation laws) from
   vector vFull to v (essentially copies indices that are not lead
   coefficients). */
void nw_cons_copyIndependentFromFull(const nw_cons *c, double *v, const double *vFull) {
  int rows=nw_cons_getNcons(c), cols=nw_cons_getNreact(c);
  int i, j, k=0;
  int lead=nw_cons_getLead(c, 0), nextLead;
  for(i=1; i<=rows; i++) {
    /*Copy everything between lead coefficient, last "lead coefficient" is the last column*/
    if(i==rows) nextLead=cols;
    else nextLead=nw_cons_getLead(c, i);
    for(j=lead+1; j<nextLead; j++) {
      /* fprintf(OUT, "copyIndependent(): Copying w[%i]=%f to v[%i], lead=%i, nextLead=%i\n",  */
      /* 	      j, w[j], k, lead, nextLead); */
      v[k++]=vFull[j];
    } 
    lead=nextLead;
  }
}

/* Copy 'independent quantities' (according to conservation laws) from
   vector vFull to v (essentially copies indices that are not lead
   coefficients). */
void nw_cons_copyIndependentToFull(const nw_cons *c, double *vFull, const double *v) {
  int rows=nw_cons_getNcons(c), cols=nw_cons_getNreact(c);
  int i, j, k=0;
  int lead=nw_cons_getLead(c, 0), nextLead;
  for(i=1; i<=rows; i++) {
    /*Copy everything between lead coefficient, last "lead coefficient" is the last column*/
    if(i==rows) nextLead=cols;
    else nextLead=nw_cons_getLead(c, i);
    for(j=lead+1; j<nextLead; j++) {
      /* fprintf(OUT, "copyIndependentToFull(): Copying w[%i]=%f to vFull[%i], lead=%i, nextLead=%i\n", */
      /* 	      k, v[k], j, lead, nextLead); */
      vFull[j]=v[k++];
    } 
    lead=nextLead;
  }
}

/* Returns derivative of DEPENDENT quantity i by independent quantity j */
double nw_cons_dX_Li_d_Xj(const nw_cons *c, int i, int j) {
  /* int indI=nw_cons_getLead(c, i); */
  /* Negative sign after solving for dependent quantity */
  /* fprintf(ERR, "nw_cons_dX_Li_d_Xj(): i=%i, indI=%i, j=%i\n", i, indI,j); */
  return -nw_cons_getKij(c,i,j);
} 



#ifdef TEST
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define NCOMP 3
#define NREACT 1

static nw_cons *consABC;

/* Calculate reduced system of par with conserved quantities removed. It is assumed that y and dydt have allocated sufficient */
int nw_cons_reducedODEs(double t,
			const double y[],
			double dydt[], void *par) {
  int err;

  /* Initialise vector of full system with independent quantities. */
  /* for(i=0; i<nC-nCons; i++) { */
  /*   yFull[i+nCons]=y[i]; */
  /* } */
  nw_cons_copyIndependentToFull(consABC, yFull, y);
  /* copyIndependentToFull(&Gview.matrix, yFull, y); */
  /* fprintf(OUT, "nw_cons_reducedODES: Copied independent quantities...\n"); */

  nw_cons_xCons(consABC, yFull);
  /* fprintf(OUT, "nw_cons_reducedODES: Calculated dependent quantities...\n"); */

  err=nw_massaction_system_ODEs(t, yFull, dydtFull, par);
  /* fprintf(OUT, "nw_cons_reducedODES: Solved full system...\n"); */


  if(err != GSL_SUCCESS) return err;

  /* for(i=0; i<nC-nCons; i++) { */
  /*   dydt[i]=dydtFull[i+nCons]; */
  /* } */
  nw_cons_copyIndependentFromFull(consABC, dydt, dydtFull);
  /* fprintf(OUT, "nw_cons_reducedODES: Copied independent derivatives...\n"); */
    
  return GSL_SUCCESS;
}

int main(int argc, char **argv) {

  /* nw_massaction *SEtoC=nw_massaction_alloc(); */
  /* nw_massaction *CtoPE=nw_massaction_alloc(); */
  double k1F=0.1, k1B=0.2;
  const char *names[]={"A", "B", "C"};
  nw_massaction_system *  mm=nw_massaction_system_alloc(NCOMP, NREACT, names);

  gsl_odeiv2_system gslModel = {nw_massaction_system_ODEs, NULL, NCOMP, mm};
  gsl_odeiv2_system  gslReducedModel = {/* nw_cons_reducedODEs */nw_massaction_system_reducedODEs, NULL, 1, mm};
  gsl_odeiv2_driver * d, *dRed;

  int i, j;
  double t = 0.0, t1 = 100.0;
  double tRed = 0.0;
  double y[NCOMP] = { 1.0, 0.5, 0.1 };
  double yRed[]={0.1};

  double Garray[2*NCOMP]={1,0,1,
			  0,1,1};

  gsl_matrix *Smatrix;

  fprintf(ERR, "Setting up reactions...\n");

  nw_massaction_system_addReact(mm, 0, 1, 2, k1F, k1B);

  fprintf(ERR, "Done...\n");

  nw_massaction_system_print(ERR, mm);

  fprintf(ERR, "Calculating stoichiometric matrix S\n");
  Smatrix=nw_massaction_system_S(mm);
  
  /* gsl_matrix_fprintf (ERR, Smatrix, "%2g"); */
  for(i=0; i<Smatrix->size1; i++) {
    for(j=0; j<Smatrix->size2; j++) {
      fprintf(ERR, "%2g\t", gsl_matrix_get(Smatrix, i, j));
    }
    fprintf(ERR, "\n");
  }

  /* consABC=nw_cons_init(Garray, 2, NCOMP); */
  /* nw_cons_conservedQuantities(consABC, y); */

  /* nw_cons_allocODE(consABC); */
  nw_massaction_system_addCons(mm, Garray, 2, NCOMP);
  nw_massaction_system_conservedQuantities(mm, y);

  d = gsl_odeiv2_driver_alloc_y_new (&gslModel, gsl_odeiv2_step_rk8pd,
				     1e-6, 1e-6, 0.0);
  dRed = gsl_odeiv2_driver_alloc_y_new (&gslReducedModel, 
					gsl_odeiv2_step_rk8pd,
				     1e-6, 1e-6, 0.0);

  /* setCons(0,1.69658); */

  /* for(i=0; i<NCOMP; i++) { */
  /*   yFull[i]=y[i]; */
  /* } */

  printf ("%g\t%g\t%g\t%g", t, y[0], y[1], y[2]);
  printf("\t%g\n", yRed[0]);

  for (i = 1; i <= 100; i++)
    {
      double ti = i * t1 / 100.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
     	{
     	  printf ("error, return value=%d\n", status);
     	  break;
     	}

      status = gsl_odeiv2_driver_apply (dRed, &tRed, ti, yRed);

      printf ("%g\t%g\t%g\t%g", t, y[0], y[1], y[2]);
      printf("\t%g\n", yRed[0]);
    }

  nw_cons_free(consABC);

  gsl_odeiv2_driver_free (d);
  gsl_odeiv2_driver_free (dRed);

  return 0;  
}
#endif
 
