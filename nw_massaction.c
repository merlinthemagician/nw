/*
 *  nw_massaction.c
 *
 * Chemical reaction network based upon mass action and first-order
 * kinetics.
 * 
 * Ivo Siekmann, 06/03/2014
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

#include <gsl/gsl_matrix.h>

#include "nw_massaction.h"

#define OUT stdout
#define ERR stderr

/* Generate nw_massaction */
nw_massaction * nw_massaction_alloc(size_t nA, size_t nB, double kF, double kB) {
  nw_massaction *r=malloc(sizeof(*r));
  r->nA=nA;
  r->A=malloc(nA*sizeof(*(r->A)));

  r->nB=nB;
  r->B=malloc(nB*sizeof(*(r->B)));

  r->kF=kF;
  r->kB=kB;

  return r;
}

/* Initialise nw_massaction for a reaction A + B -> C*/
nw_massaction *nw_massaction_ABC(int i, int j, int k, 
				 char *nameA, char *nameB, char* nameC,
				 double kF, double kB) {
  nw_massaction *r=nw_massaction_alloc(2,1,kF,kB);

  nw_massaction_setA(r, 0, nameA, i, 1);
  nw_massaction_setA(r, 1, nameB, j, 1);

  nw_massaction_setB(r, 0, nameC, k, 1);

  return r;
}

/* Initialise nw_massaction for a reaction A + B -> A + C*/
/* No backward reaction */
nw_massaction *nw_massaction_enzyme(int i, int j, int k, char *nameA, char *nameB, char* nameC,
				    double kF) {
  nw_massaction *r=nw_massaction_alloc(2,2,kF,0);

  nw_massaction_setA(r, 0, nameA, i, 1);
  nw_massaction_setA(r, 1, nameB, j, 1);

  nw_massaction_setB(r, 0, nameA, i, 1);
  nw_massaction_setB(r, 1, nameC, k, 1);

  return r;
}

/* Get reactant i from complex A */
nw_reactant *nw_massaction_getAReactant(const nw_massaction *r, int i) {
  return r->A[i];
}

/* Get index of reactant i from complex A */
int nw_massaction_getAr(const nw_massaction *r, int i) {
  nw_reactant *react=r->A[i];

  return nw_reactant_getR(react);
}

/* Get stoichmetric coefficient of reactant i from complex A */
int nw_massaction_getAs(const nw_massaction *r, int i) {
  nw_reactant *react=r->A[i];

  return nw_reactant_getStoichiometry(react);
}

/* Set reactant i in complex A */
void nw_massaction_setAReactant(const nw_massaction *r, int i, nw_reactant *react) {
  r->A[i]=react;
}

/* Set reactant i in complex A */
void nw_massaction_setA(const nw_massaction *r, int i, 
			const char *name, int Ak, int As) {
  nw_reactant *react=nw_reactant_alloc(name, Ak, As);
  nw_massaction_setAReactant(r, i, react);
}


/* Get reactant i from complex B */
nw_reactant *nw_massaction_getBReactant(const nw_massaction *r, int i) {
  return r->B[i];
}

/* Set reactant i in complex B */
void nw_massaction_setBReactant(const nw_massaction *r, int i, nw_reactant *react) {
  r->B[i]=react;
}

/* Set reactant i in complex B */
void nw_massaction_setB(const nw_massaction *r, int i, 
			const char *name, int Bk, int Bs) {
  nw_reactant *react=nw_reactant_alloc(name, Bk, Bs);
  nw_massaction_setBReactant(r, i, react);
}

/* Get index of reactant i from complex A */
int nw_massaction_getBr(const nw_massaction *r, int i) {
  nw_reactant *react=r->B[i];

  return nw_reactant_getR(react);
}

/* Get stoichmetric coefficient of reactant i from complex B */
int nw_massaction_getBs(const nw_massaction *r, int i) {
  nw_reactant *react=r->B[i];

  return nw_reactant_getStoichiometry(react);
}

/* Returns forward rate kF. */
double nw_massaction_get_kF(const nw_massaction *r) {
  return r->kF;
}

/* Returns backwards rate kB. */
double nw_massaction_get_kB(const nw_massaction *r) {
  return r->kB;
}

/* Sets rate forward rate kF. */
void nw_massaction_set_kF(nw_massaction *r,  double kF) {
  r->kF=kF;
}

/* Sets backward rate kB. */
void nw_massaction_set_kB(nw_massaction *r, double kB) {
  r->kB=kB;
}

void nw_massaction_printComplex(FILE *fp, const nw_reactant **C, size_t nC) {
  char *strCi;
  int i;
  for(i=0; i<nC-1; i++) {
    strCi=nw_reactant_toString(C[i]);
    fprintf(fp, "%s + ", strCi);
    /* free(strCi); */
  }
  strCi=nw_reactant_toString(C[i]);
  fprintf(fp, "%s ", strCi);
  /* free(strCi); */
}

/* Print phosphorylation network */
void nw_massaction_print(FILE *fp, const nw_massaction *r) {
  nw_massaction_printComplex(fp, (const nw_reactant **)r->A, r->nA);
  fprintf(fp, " <-> ");
  nw_massaction_printComplex(fp, (const nw_reactant **)r->B, r->nB);

  fprintf(fp, ",\tkF=%g, kB=%g\n", 
	  nw_massaction_get_kF(r),
	  nw_massaction_get_kB(r));
}

/* Mass action law for complex C containing nC reactants with concentrations y */
double nw_massaction_ComplexReaction(const nw_reactant **C, size_t nC, 
				     const double *y) {
  double prod=1;
  int i;

  for(i=0; i<nC; i++) {
    int k=nw_reactant_getR(C[i]), s=nw_reactant_getStoichiometry(C[i]);
    double yS=pow(y[k], s);
    /* if (yS==0) return 0; */
    prod *=yS;
  }
  return prod;
}

/* Mass action reaction */
double nw_massaction_reaction(const nw_massaction *r, const double *y) {
  nw_reactant **A=r->A, **B=r->B;
  double compA, compB;
  double kF=r->kF, kB=r->kB;

  compA=nw_massaction_ComplexReaction((const nw_reactant **)A, r->nA, y);
  compB=nw_massaction_ComplexReaction((const nw_reactant **)B, r->nB, y);

  return -kF*compA + kB*compB;
}

/* Add contribution of reaction r to dydt */
void nw_massaction_setODE(const nw_massaction *r, const double *y, double *dydt) {
  int i;
  double react=nw_massaction_reaction(r, y);

  /* left-hand side */
  for(i=0; i<r->nA; i++) {
    int k=nw_massaction_getAr(r, i);
    dydt[k] += react;
  }

  /* right-hand side */
  for(i=0; i<r->nB; i++) {
    int k=nw_massaction_getBr(r, i);
    dydt[k] -= react;
  }
}


/* Partial derivative by x_n of complex C containing nC reactants with
   concentrations y */
double nw_massaction_ComplexReactionDxN(const nw_reactant **C, size_t nC, 
					const double *y, int n) {
  double prod=1;
  int i;
  /* a priori, derivative is zero */
  int isZero=1;

  /* Quite inefficient because it is not clear until the end if x_n is part of the product! */
  for(i=0; i<nC; i++) {
    int k=nw_reactant_getR(C[i]), s=nw_reactant_getStoichiometry(C[i]);
    double yS;
    if(k!=n){
      yS=pow(y[k], s);
    }
    else {
      yS=s*pow(y[k], s-1), isZero=0;
    }
    /* if (yS==0) return 0; */
    prod *=yS;
  }
  if(isZero) return 0;

  return prod;
}

/* Derivative of Mass action reaction */
double nw_massaction_dxN(const nw_massaction *r, const double *y, int n) {
  nw_reactant **A=r->A, **B=r->B;
  double DcompADxn, DcompBDxn;
  double kF=r->kF, kB=r->kB;

  DcompADxn=nw_massaction_ComplexReactionDxN((const nw_reactant **)A, r->nA, y, n);
  DcompBDxn=nw_massaction_ComplexReactionDxN((const nw_reactant **)B, r->nB, y, n);

  return -kF*DcompADxn + kB*DcompBDxn;
}

/* Add contribution of reaction r to dfdy */
void nw_massaction_setJac(const nw_massaction *r, const double *y, double *dfdy, int nC) {
  int i, n;
  gsl_matrix_view dfdy_mat= gsl_matrix_view_array (dfdy, nC, nC);
  gsl_matrix * J = &dfdy_mat.matrix;

  /* left-hand side... */
  for(i=0; i<r->nA; i++) {
    int k=nw_massaction_getAr(r, i);
    for(n=0; n<nC; n++) {
      double Jkn=gsl_matrix_get(J, k, n);
      /* ... positive contribution */
      Jkn+=nw_massaction_dxN(r, y, n);
      gsl_matrix_set(J, k, n, Jkn);
    }
  }

  /* right-hand side... */
  for(i=0; i<r->nB; i++) {
    int k=nw_massaction_getBr(r, i);
    for(n=0; n<nC; n++) {
      double Jkn=gsl_matrix_get(J, k, n);
      /* ... negative contribution */
      Jkn-=nw_massaction_dxN(r, y, n);
      gsl_matrix_set(J, k, n, Jkn);
    }
  }
}

/* Set components of Jacobian for mass action reaction r. It
   is assumed that the system has nR reactions and nC reactants. The
   array dfdy contains matrix elements in column-first order. */
/* void nw_massaction_setJacR(const nw_massaction *r, const double *y,  */
/* 			   double *dfdy, int nC) { */
/*   gsl_matrix_view dfdy_mat= gsl_matrix_view_array (dfdy, nC, nC); */
/*   gsl_matrix * J = &dfdy_mat.matrix; */
/*   int k, n; */

/*   /\* Differentiate reaction r by all possible n *\/ */
/*   for(k=0; k<nC; k++) { */
/*     for(n=0; n<nC; n++) { */
/*       double Jkn=gsl_matrix_get(J, k, n); */
/*       Jkn+=nw_massaction_dxN(r, y, n); */
/*       gsl_matrix_set(J, k, n, Jkn); */
/*     } */
/*   } */
/* } */

/* Set components of Jacobian for the k-th mass action reaction r. It
   is assumed that the system has nR reactions and nC reactants. The
   array dfdy contains matrix elements in column-first order. */
/* void nw_massaction_setJac(const nw_massaction *r, const double *y,  */
/* 			    double *dfdy, int nR, int nC) { */
/*   gsl_matrix_view dfdy_mat= gsl_matrix_view_array (dfdy, nR, nC); */
/*   gsl_matrix * J = &dfdy_mat.matrix; */
/*   int k; */

/*   for(k=0; k<nR; k++) { */
/*     nw_massaction_setJacRk(r, y, dfdy, k, nR, nC); */
/*   } */
/* } */
