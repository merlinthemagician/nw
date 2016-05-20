#ifndef NW_MASSACTION
#define NW_MASSACTION
/*
 *  nw_massaction.h
 *
 * Chemical reaction network based upon mass action kinetics
 * 
 * Ivo Siekmann, 07/03/2014
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include "nw_reactant.h"

/* A and B are complexes of reactants. 
   kF, kB: Forward and backward rate. */
typedef struct {
  /*Complex A*/
  nw_reactant **A;
  size_t nA;

  /*Complex B*/
  nw_reactant **B;
  size_t nB;

  /* forward rate */
  double kF;
  /* backward rate */
  double kB;
} nw_massaction;

/*Generate nw_massaction*/
nw_massaction * nw_massaction_alloc(size_t nA, size_t nB, double kF, double kB);

/* Sets mass action reaction at (i,j) and labels product.*/
/* void nw_massaction_setParams(nw_model *m, int i, int j, nw_massaction *e,  */
/* 			     label product); */

/* Get parameters for network model at index k */
/* nw_massaction* nw_massaction_getParams(const nw_model *m, int k); */

/* Initialise nw_massaction for a reaction A + B -> C*/
nw_massaction *nw_massaction_ABC(int i, int j, int k, 
				 char *A, char *B, char* C,
				 double kF, double kB);

/* Initialise nw_massaction for a reaction A + B -> A + C*/
/* No backward reaction */
nw_massaction *nw_massaction_enzyme(int i, int j, int k, 
				    char *A, char *B, char* C,
				    double kF);

/* Get reactant i from complex A */
nw_reactant *nw_massaction_getAReactant(const nw_massaction *r, int i);

/* Get index of reactant i from complex A */
int nw_massaction_getAr(const nw_massaction *r, int i);

/* Get stoichometric coefficient of reactant i from complex A */
int nw_massaction_getAs(const nw_massaction *r, int i);

/* Set reactant i in complex A */
void nw_massaction_setAReactant(const nw_massaction *r, int i, nw_reactant *react);

/* Set reactant i in complex A */
void nw_massaction_setA(const nw_massaction *r, int i, const char *name, int Ak, int As);

/* Get reactant i from complex B */
nw_reactant *nw_massaction_getBReactant(const nw_massaction *r, int i);

/* Get index of reactant i from complex B */
int nw_massaction_getBr(const nw_massaction *r, int i);

/* Get stoichometric coefficient of reactant i from complex B */
int nw_massaction_getBs(const nw_massaction *r, int i);

/* Set reactant i in complex B */
void nw_massaction_setBReactant(const nw_massaction *r, int i, nw_reactant *react);

/* Set reactant i in complex B */
void nw_massaction_setB(const nw_massaction *r, int i, const char *name, int Bk, int Bs);

/* Returns forward rate kF. */
double nw_massaction_get_kF(const nw_massaction *r);

/* Returns backwards rate kB. */
double nw_massaction_get_kB(const nw_massaction *r);

/* Returns pointer to forward rate kF_ij. */
/* double * nw_massaction_getModel_kF_Ptr(const nw_model *m, int i, int j); */

/* Returns forward rate kF_ij for protein k. */
/* double nw_massaction_getModel_kF(const nw_model *m, int i, int j); */

/* Returns backward rate kB_ij. */
/* double * nw_massaction_getModel_kB_Ptr(const nw_model *m, int k, int i); */

/* Sets rate forward rate kF. */
void nw_massaction_set_kF(nw_massaction *r,  double kF);

/* Sets forward rate kF_ij. */
/* void nw_massaction_setModel_kF(nw_model *m, int k, int i, double Vi); */

/* Sets backward rate kB_ij. */
void nw_massaction_set_kB(nw_massaction *r, double kB);

/* Print phosphorylation network */
void nw_massaction_print(FILE *fp, const nw_massaction *r);

/* Add contribution of reaction r to dydt */
void nw_massaction_setODE(const nw_massaction *r, const double *y, double *dydt);

/* Set components of Jacobian for mass action reaction r. It is
   assumed that the system has nR reactions and nC reactants. The
   array dfdy contains matrix elements in column-first order. */
void nw_massaction_setJac(const nw_massaction *r, const double *y, 
			  double *dfdy, int nC);
#endif
