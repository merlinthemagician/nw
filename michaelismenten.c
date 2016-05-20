/*
 * michaelismenten.c
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

#define NCOMP 4
#define NREACT 2

#define OUT stdout
#define ERR stderr

static double Garray[]={1,0,1,1,
			0,1,1,0};

static int enzymeApp(double t, const double y[], double dydt[], void *par) {
  double* p=(double *)par;
  dydt[0] = -p[0]*y[0]/(y[0] + p[1]);

  return GSL_SUCCESS;
}

int main(int argc, char **argv) {

  /* nw_massaction *SEtoC=nw_massaction_alloc(); */
  /* nw_massaction *CtoPE=nw_massaction_alloc(); */
  double k1F=100, k1B=0.1;
  double k2F=0, k2B=1.5;
  const char *names[]={"S", "E", "C", "P"};
  nw_massaction_system *  mm=nw_massaction_system_alloc(NCOMP, NREACT, names);

  gsl_odeiv2_system gslModel = {nw_massaction_system_ODEs, NULL, NCOMP, mm};
  gsl_odeiv2_driver * d;

  double par[2];
  gsl_odeiv2_system gslEnzyme = {enzymeApp, NULL, 1, par};
  gsl_odeiv2_driver * d1;

  int i, j;
  double t = 0.0, t1 = 100.0;
  double tau=0.0;
  double y[NCOMP] = { 1.0, 0.01, 0.0, 0.0 };
  double S=y[0];

  gsl_matrix *Smatrix;

  double *Jarray=malloc(NCOMP*NCOMP*sizeof(*Jarray));
  gsl_matrix_view J=gsl_matrix_view_array(Jarray, NCOMP, NCOMP);

  double *Jredarray=malloc((NCOMP-2)*(NCOMP-2)*sizeof(*Jredarray));
  gsl_matrix_view Jred=gsl_matrix_view_array(Jredarray, NCOMP-2, NCOMP-2);

  double yRed[NCOMP-2];

  /* nw_massaction_init(SEtoC, 0, 1, 2, k1F, k1B); */
  /* nw_massaction_init(CtoPE, 1, 3, 2, k2F, k2B); */

  fprintf(ERR, "Adding reactions...\n");

  nw_massaction_system_addReact(mm, 0, 1, 2, k1F, k1B);
  nw_massaction_system_addReact(mm, 1, 3, 2, k2F, k2B);

  fprintf(ERR, "Done...\n");

  /* mm->names[0]="S"; */
  /* mm->names[1]="E"; */
  /* mm->names[2]="C"; */
  /* mm->names[3]="P"; */

  /* nw_massaction_system_setReact(mm, 0, SEtoC); */
  /* nw_massaction_system_setReact(mm, 1, CtoPE);   */

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

  /* Parameters from approximation */
  par[0]=k2B/k1F;
  par[1]=(k1B + k2B)/k1F;

  d = gsl_odeiv2_driver_alloc_y_new (&gslModel, gsl_odeiv2_step_rk8pd,
				     1e-6, 1e-6, 0.0);

  d1 = gsl_odeiv2_driver_alloc_y_new (&gslEnzyme, gsl_odeiv2_step_rk8pd,
				      1e-6, 1e-6, 0.0);

  fprintf(stdout, "t\t");
  nw_massaction_system_printNames(stdout, mm);

  nw_massaction_system_addCons(mm, Garray, 2, NCOMP);
  nw_massaction_system_conservedQuantities(mm, y);

  printf ("%g\t%g\t%g\t%g\t%g", t, y[0], y[1], y[2], y[3]);
  printf("\t%g\n", S);

  for (i = 1; i <= 100; i++)
    {
      double ti = i * t1 / 100.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
     	{
     	  printf ("error, return value=%d\n", status);
     	  break;
     	}

      gsl_odeiv2_driver_apply (d1, &tau, ti, &S);
      
      if (status != GSL_SUCCESS)
     	{
     	  printf ("error, return value=%d\n", status);
     	  break;
     	}

      printf ("%g\t%g\t%g\t%g\t%g", t, y[0], y[1], y[2], y[3]);
      printf("\t%g\n", S);
    }

  gsl_odeiv2_driver_free (d);
  gsl_odeiv2_driver_free (d1);


  fprintf (ERR, "Jacobian:\n");
  nw_massaction_system_setJac(mm, y, Jarray);
  for(i=0; i<J.matrix.size1; i++) {
    for(j=0; j<J.matrix.size2; j++) {
      fprintf(ERR, "%g\t", gsl_matrix_get(&J.matrix, i, j));
    }
    fprintf(ERR, "\n");
  }

  /* fprintf (ERR, "Reduced Jacobian (nw_massaction_totalDerivative()):\n"); */
  /* fprintf(ERR, "J33=%g\n", nw_massaction_system_totalDerivative(mm, &J.matrix, 2, 2)); */

  yRed[0]=y[2];
  yRed[1]=y[3];

  nw_massaction_system_setReducedJac(mm, yRed, Jredarray);
  fprintf (ERR, "Reduced Jacobian (nw_massaction_system_reducedJac()):\n");
  for(i=0; i<Jred.matrix.size1; i++) {
    for(j=0; j<Jred.matrix.size2; j++) {
      fprintf(ERR, "%g\t", gsl_matrix_get(&Jred.matrix, i, j));
    }
    fprintf(ERR, "\n");
  }
  return 0;  
}
