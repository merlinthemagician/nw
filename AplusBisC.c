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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "nw_massaction.h"
#include "nw_massaction_system.h"

#define NCOMP 3
#define NREACT 1

#define OUT stdout
#define ERR stderr



#ifdef TEST
static double Garray[]={1,0,1,
			0,1,1};


int main(int argc, char **argv) {

  /* nw_massaction *SEtoC=nw_massaction_alloc(); */
  /* nw_massaction *CtoPE=nw_massaction_alloc(); */
  double k1F=5, k1B=0.2;
  const char *names[]={"A", "B", "C"};
  nw_massaction_system *  mm=nw_massaction_system_alloc(NCOMP, NREACT, names);

  gsl_odeiv2_system gslModel = {nw_massaction_system_ODEs, NULL, NCOMP, mm};
  gsl_odeiv2_system  gslReducedModel = {nw_massaction_system_reducedODEs, 
					NULL, 1, mm};
  gsl_odeiv2_driver * d, *dRed;

  int i, j;
  double t = 0.0, t1 = 100.0;
  double tRed = 0.0;
  double y[NCOMP] = { 1.0, 0.5, 0.1 };
  double yRed[]={0.1};

  double J[NCOMP*NCOMP];
  gsl_matrix_view Jview=gsl_matrix_view_array(J, NCOMP, NCOMP);

  gsl_matrix *Smatrix;

  fprintf(ERR, "Setting up reactions...\n");

  nw_massaction_system_addReact(mm, 0, 1, 2, k1F, k1B);

  fprintf(ERR, "Done...\n");

  nw_massaction_system_print(ERR, mm);

  fprintf(ERR, "Calculating stoichiometric matrix S\n");
  Smatrix=nw_massaction_system_S(mm);

  nw_massaction_system_addCons(mm, Garray, 2, NCOMP);
  nw_massaction_system_conservedQuantities(mm, y);

  
  fprintf (ERR, "Stoichiometric matrix S:\n");
  for(i=0; i<Smatrix->size1; i++) {
    for(j=0; j<Smatrix->size2; j++) {
      fprintf(ERR, "%2g\t", gsl_matrix_get(Smatrix, i, j));
    }
    fprintf(ERR, "\n");
  }

  fprintf (ERR, "Jacobian:\n");
  nw_massaction_system_setJac(mm, y, J);
  for(i=0; i<Jview.matrix.size1; i++) {
    for(j=0; j<Jview.matrix.size2; j++) {
      fprintf(ERR, "%g\t", gsl_matrix_get(&Jview.matrix, i, j));
    }
    fprintf(ERR, "\n");
  }

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

  gsl_odeiv2_driver_free (d);
  gsl_odeiv2_driver_free (dRed);

  fprintf (ERR, "Jacobian:\n");
  nw_massaction_system_setJac(mm, y, J);
  for(i=0; i<Jview.matrix.size1; i++) {
    for(j=0; j<Jview.matrix.size2; j++) {
      fprintf(ERR, "%g\t", gsl_matrix_get(&Jview.matrix, i, j));
    }
    fprintf(ERR, "\n");
  }

  fprintf (ERR, "Reduced Jacobian (nw_massaction_totalDerivative()):\n");
  fprintf(ERR, "J33=%g\n", nw_massaction_system_totalDerivative(mm, J, 2, 2));

  double J33=0;
  /* Set Jacobian for reduced reactions of massaction_system sys. */
  nw_massaction_system_setReducedJac(mm, yRed, &J33);
  fprintf (ERR, "Reduced Jacobian (nw_massaction_system_reducedJac()):\n");
  fprintf(ERR, "J33=%g\n", J33);
  /* nw_massaction_system_setJac(mm, y, J); */
  /* for(i=0; i<Jview.matrix.size1; i++) { */
  /*   for(j=0; j<Jview.matrix.size2; j++) { */
  /*     fprintf(ERR, "%g\t", gsl_matrix_get(&Jview.matrix, i, j)); */
  /*   } */
  /*   fprintf(ERR, "\n"); */
  /* } */


  return 0;  
}
#endif
