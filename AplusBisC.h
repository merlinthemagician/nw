#ifndef APLUSBISC
#define APLUSBISC
/*
 * AplusBisC.c
 *
 *
 * Simple mass action kinetics: A + B <-> C
 *
 * Ivo Siekmann, 29/04/2014
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#include <stdio.h>

#define NCOMP 3
#define NREACT 1

#define OUT stdout
#define ERR stderr

/* int nw_massaction_system_ODEs(double t, const double y[], double dydt[], void *par); */

/* Set conserved quantity i to C */
/* void setCons(int i, double C); */

/* Get conserved quantity i */
/* double getCons(int i); */

/* Calculate concentration of n-th reactant from conserved quantity T */
/* double xConsN(const gsl_matrix *G, const gsl_vector *T, */
/* 	      const gsl_vector* x, int n); */

/*Calculate all conserved quantities and save to first nC components
  of y (it should be safe to do this if G is in echelon form?!?) */
/* void xCons(const gsl_matrix *G, const gsl_vector *T, */
/* 	   gsl_vector* x); */

/* Calculate reduced system with first nCons conserved quantities
   removed */
/* int nw_massaction_system_reducedODEs(double t,  */
/* 				     const double y[],  */
/* 				     double dydt[], void *par); */
#endif
