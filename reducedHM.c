/*
 *  HockinMann_JBC2002.c
 *
 *
 * Hockin et al. & Mann coagulation model
 *
 *
 * Ivo Siekmann, 20/03/2014
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
#include "nw_cons.h"

#include "HockinMann_JBC2002.h"

#define OUT stdout
#define ERR stderr

int main(int argc, char **argv) {

  nw_massaction_system *  hm_jbc2002=init_hm_jbc2002();

  gsl_odeiv2_system gslModel = {nw_massaction_system_ODEs/* hm2002 */,
				nw_massaction_system_jac/* jac_hm2002jac */, NCOMP, hm_jbc2002};

  gsl_odeiv2_system gslReducedModel = {nw_massaction_system_reducedODEs/* hm2002 */,
				       NULL, NCOMP, hm_jbc2002};
  gsl_odeiv2_driver * d, *dRed;
  double y[NCOMP] = { 5*1e-3, 
		      VII0, 
		      0, 
		      VIIa0, 
		      0,
		      0,
		      0,
		      X0, 
		      0,
		      0,
		      IX0,
		      0,
		      0,
		      II0,
		      VIII0,
		      0,
		      0,
		      0,
		      0,
		      0,
		      V0, 
		      0,
		      0,
		      0,
		      0,
		      TFPI0,
		      0,
		      0,
		      ATIII0,
		      0,
		      0,
		      0,
		      0,
		      0
};

  double yRed[NCOMP-10];

  int i, j;
  double t=0, tRed=0, t1=100, tend=700/* 22279 *//* 10000 */;
  int nBlood=34;
  gsl_matrix *Smatrix;
  int withTFPI=0, withATIII=0;
  double *conserved;

  double J[NCOMP*NCOMP];
  gsl_matrix_view Jview=gsl_matrix_view_array(J, NCOMP, NCOMP);

  double Jred[(NCOMP-2)*(NCOMP-2)];
  gsl_matrix_view Jredview=gsl_matrix_view_array(Jred, NCOMP-2, NCOMP-2);

  nw_massaction_system_print(ERR, hm_jbc2002);

  /*   fprintf(ERR, "Calculating stoichiometric matrix S\n"); */
  /* Smatrix=nw_massaction_system_S(hm_jbc2002); */
  
  /* gsl_matrix_fprintf (ERR, Smatrix, "%2g"); */
  /* for(i=0; i<Smatrix->size1; i++) { */
  /*   for(j=0; j<Smatrix->size2; j++) { */
  /*     fprintf(OUT, "%2g\t", gsl_matrix_get(Smatrix, i, j)); */
  /*   } */
  /*   fprintf(OUT, "\n"); */
  /* } */

  /* exit(1); */
  nw_massaction_system_addCons(hm_jbc2002, kerHM, 10, nBlood);

  /*ARGUMENTS*/
  if(argc>=2) {
    y[0]=atof(argv[1])*1e-3;
  }

  if(argc>=3) {
    withTFPI=atoi(argv[2]);
  }

  if(argc>=4) {
    withATIII=atoi(argv[3]);
  }
  
  /* Pro-Coagulation model only - set TFPI and ATPIII to zero!*/
  if(! withTFPI) y[25]=0;
  if(! withATIII) y[28]=0;

  conserved=nw_massaction_system_conservedQuantities(hm_jbc2002, y);
  fprintf(OUT, "Conserved quantities:\n");
  for(i=0; i<10; i++) {
    fprintf(OUT, "T_%i=%g\n", i, conserved[i]);
  }

  nw_massaction_system_copyIndependentFromFull(hm_jbc2002, yRed, y);

  d = gsl_odeiv2_driver_alloc_y_new (&gslModel, gsl_odeiv2_step_rkf45/* gsl_odeiv2_step_msbdf *//* gsl_odeiv2_step_bsimp *//* gsl_odeiv2_step_rk8pd */,
				     1e-6, 1e-6, 0.0);
  dRed = gsl_odeiv2_driver_alloc_y_new (&gslReducedModel, gsl_odeiv2_step_rkf45/* gsl_odeiv2_step_msbdf *//* gsl_odeiv2_step_bsimp *//* gsl_odeiv2_step_rk8pd */,
				     1e-6, 1e-6, 0.0);

  /* Title */
  fprintf(OUT, "t\t");
  /* for(i=0; i<nBlood-1; i++) fprintf(OUT, "%s\t", nw_massaction_system_getName(hm_jbc2002, i)); */
  /* fprintf(OUT, "%s\t", nw_massaction_system_getName(hm_jbc2002, i)); */
  for(i=0; i<nBlood; i++) fprintf(OUT, "%s\t", nw_massaction_system_getName(hm_jbc2002, i));
  /* fprintf(OUT, "%s\t", nw_massaction_system_getName(hm_jbc2002, i)); */
/* Print names of reactants in massaction_system sys */
  /* nw_massaction_system_printNames(OUT, hm_jbc2002); */
  /* fprintf(OUT, "\t"); */
  nw_massaction_system_printIndependentNames(OUT, hm_jbc2002);

  fprintf(OUT, "%g\t", 0.0);
  for(i=0; i<nBlood-1; i++) fprintf(OUT, "%g\t", y[i]);
  fprintf(OUT, "%g\t", y[i]);

  for(i=0; i<nBlood-10-1; i++) fprintf(OUT, "%g\t", yRed[i]);
  fprintf(OUT, "%g\n", yRed[i]);

  /* Solving the model */
  for (i = 1; i <= tend; i++)
    {
      double ti = i /* * tend / 100.0 */;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
      int k;

      if (status != GSL_SUCCESS)
     	{
     	  printf ("error, return value=%d\n", status);
     	  break;
     	}

      fprintf(OUT, "%g\t", t);
      for(k=0; k<nBlood-1; k++) fprintf(OUT, "%g\t", y[k]);
      fprintf(OUT, "%g\n", y[k]);

      status = gsl_odeiv2_driver_apply (dRed, &tRed, ti, yRed);

      if (status != GSL_SUCCESS)
     	{
     	  printf ("error, return value=%d\n", status);
     	  break;
     	}
      for(k=0; k<nBlood-10-1; k++) fprintf(OUT, "%g\t", yRed[k]);
      fprintf(OUT, "%g\n", yRed[k]);
    }


  gsl_odeiv2_driver_free (d);
  gsl_odeiv2_driver_free (dRed);


  /* fprintf (ERR, "Jacobian:\n"); */
  /* nw_massaction_system_setJac(hm_jbc2002, y, J); */
  /* for(i=0; i<Jview.matrix.size1; i++) { */
  /*   for(j=0; j<Jview.matrix.size2; j++) { */
  /*     fprintf(ERR, "%g\t", gsl_matrix_get(&Jview.matrix, i, j)); */
  /*   } */
  /*   fprintf(ERR, "\n"); */
  /* } */

  /* nw_massaction_system_reducedJac(hm_jbc2002, yRed, Jred); */
  /* fprintf (ERR, "Reduced Jacobian (nw_massaction_system_reducedJac()):\n"); */
  /* for(i=0; i<Jredview.matrix.size1; i++) { */
  /*   for(j=0; j<Jredview.matrix.size2; j++) { */
  /*     double Jij=gsl_matrix_get(&Jredview.matrix, i, j); */
  /*     if(fabs(Jij)<1e-12) Jij=0; */
  /*     fprintf(ERR, "%g\t", Jij); */
  /*   } */
  /*   fprintf(ERR, "\n"); */
  /* } */
  

  return 0;  
}
