/*
 *  JonesMann_JBC1994.c
 *
 *
 * Jones & Mann coagulation model
 *
 * Set up nw_massaction structures and solve system of ODEs with GSL 
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

#include "nw_massaction.h"
#include "nw_massaction_system.h"

#include "JonesMann_JBC1994.h"

#define OUT stdout
#define ERR stderr

static int jm1994(double t, const double y[], double dydt[], void *par) {
  int err=nw_massaction_system_ODEs(t, y, dydt, par);
  double I=y[18];
  double VIIIaIXa=y[6];

  if(err != GSL_SUCCESS) return err;

  dydt[18]+=K20*(I-VIIIaIXa)-K20*fabs(I-VIIIaIXa);
  dydt[6]+=-fabs(I-VIIIaIXa) + (I-VIIIaIXa);

  return GSL_SUCCESS;
}


int main(int argc, char **argv) {

  nw_massaction_system *  jm_jbc1994=nw_massaction_system_alloc(NCOMP, NREACT, compNames);

  gsl_odeiv2_system gslModel = {jm1994, NULL, NCOMP, jm_jbc1994};
  gsl_odeiv2_driver * d;
  double y[NCOMP] = { 0.5, IX0, X0, V0, VIII0, II0, 0, 0, 0, 0,
		      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		      0};

  int i;
  double t=0, t1=100;
  int nBlood=19;

  /* nw_massaction_system_setNames(jm_jbc1994, compNames,NCOMP); */

  /*0*/
  nw_massaction_system_addReact(jm_jbc1994, 1, 0, 11, K6, K16);
  /*1*/
  nw_massaction_system_addReact(jm_jbc1994, 0, 14, 11, 0, K11);

  /*2*/
  nw_massaction_system_addReact(jm_jbc1994, 2, 0, 12, K6, K17);
  /*3*/
  nw_massaction_system_addReact(jm_jbc1994, 0, 15, 12, 0, K12);

  /*4*/
  nw_massaction_system_addReact(jm_jbc1994, 2, 6, 13, K6, K18);
  /*5*/
  nw_massaction_system_addReact(jm_jbc1994, 6, 15, 13, 0, K13);

  /*6*//* One reaction written down as two...*/
  nw_massaction_system_addReact(jm_jbc1994, 1, 15, 15, K15, 0);
  nw_massaction_system_addReact(jm_jbc1994, 1, 15, 14, K15, 0);

  /*7*//* One reaction written down as two...*/
  nw_massaction_system_addReact(jm_jbc1994, 3, 15, 15, K1, 0);
  nw_massaction_system_addReact(jm_jbc1994, 3, 15, 16, K1, 0);

  /*8*//* One reaction written down as two...*/
  nw_massaction_system_addReact(jm_jbc1994, 4, 15, 15, K3, 0);
  nw_massaction_system_addReact(jm_jbc1994, 4, 15, 17, K3, 0);

  /*9*//* One reaction written down as two...*/
  nw_massaction_system_addReact(jm_jbc1994, 3, 8, 8, K2, 0);
  nw_massaction_system_addReact(jm_jbc1994, 3, 8, 16, K2, 0);

  /*10*//* One reaction written down as two...*/
  nw_massaction_system_addReact(jm_jbc1994, 4, 8, 8, K4, 0);
  nw_massaction_system_addReact(jm_jbc1994, 4, 8, 17, K4, 0);

  /*11*/
  nw_massaction_system_addReact(jm_jbc1994, 5, 7, 9, K6, K19);
  nw_massaction_system_addReact(jm_jbc1994, 7, 10, 9, 0, K14);

  /*12*/
  /* One reaction written down as two...*/
  nw_massaction_system_addReact(jm_jbc1994, 10, 7, 7, K5, 0);
  nw_massaction_system_addReact(jm_jbc1994, 10, 7, 8, K5, 0);

  /*13*/
  nw_massaction_system_addReact(jm_jbc1994, 17, 14, 6, K7, K9);

  /*14*//*Two complexes must decay in order to form two individual activated comunds*/
  nw_massaction_system_addReact(jm_jbc1994, 16, 15, 7, K8, 2*K10);

  nw_massaction_system_print(ERR, jm_jbc1994);

  d = gsl_odeiv2_driver_alloc_y_new (&gslModel, gsl_odeiv2_step_rk8pd,
				     1e-6, 1e-6, 0.0);

  /* Title */
  fprintf(OUT, "t\t");
  for(i=0; i<nBlood-1; i++) fprintf(OUT, "%s\t", nw_massaction_system_getName(jm_jbc1994, i));
  fprintf(OUT, "%s\t", nw_massaction_system_getName(jm_jbc1994, i));

  fprintf(OUT, "%g\t", 0.0);
  for(i=0; i<nBlood-1; i++) fprintf(OUT, "%g\t", y[i]);
  fprintf(OUT, "%g\n", y[i]);

  /* Solving the model */
  for (i = 1; i <= 250; i++)
    {
      double ti = i * t1 / 100.0;
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
    }

  gsl_odeiv2_driver_free (d);
  

  return 0;  
}
