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

#include "HockinMann_JBC2002.h"

#define OUT stdout
#define ERR stderr

/* Initialise Hockin/Mann system */
nw_massaction_system* init_hm_jbc2002() {
  nw_massaction_system *  hm_jbc2002=nw_massaction_system_alloc(NCOMP, NREACT,compNames);

  /* nw_massaction_system_setNames(hm_jbc2002, compNames,NCOMP); */

  /*1*/ /*TF + VII <-> TF.VII */
  nw_massaction_system_addReact(hm_jbc2002, 0, 1, 2, K2, K1);

  /*2*/ /*TF + VIIa <-> TF.VIIa */
  nw_massaction_system_addReact(hm_jbc2002, 0, 3, 4, K4, K3);

  /*3*/ /*TF.VIIa + VII <-> TF.VIIa + VIIa*/
  /* nw_massaction_system_addReact(hm_jbc2002, 4, 1, 3, K5, 0); */
  /* nw_massaction_system_addReact(hm_jbc2002, 4, 1, 4, K5, 0); */
  nw_massaction_system_addEnzyme(hm_jbc2002, 4, 1, 3, K5);

  /*4*/ /*Xa + VII -> Xa + VIIa*/
  /* nw_massaction_system_addReact(hm_jbc2002, 5, 1, 3, K6, 0); */
  /* nw_massaction_system_addReact(hm_jbc2002, 5, 1, 5, K6, 0); */
  nw_massaction_system_addEnzyme(hm_jbc2002, 5, 1, 3, K6);

  /*5*/ /*IIa + VII -> IIa + VIIa */
  /* nw_massaction_system_addReact(hm_jbc2002, 6, 1, 3, K7, 0); */
  /* nw_massaction_system_addReact(hm_jbc2002, 6, 1, 6, K7, 0); */
  nw_massaction_system_addEnzyme(hm_jbc2002, 6, 1, 3, K7);


  /*6*/ /* TF.VIIa + X <-> TF.VIIa.X*/
  nw_massaction_system_addReact(hm_jbc2002, 4, 7, 8, K9, K8);
  /* FIRST - ORDER reaction */
  /* nw_massaction_system_addReact(hm_jbc2002, 6, 1, 3, K6, 0); */

  /*7*/ /*TF.VIIa + Xa <-> TF.VIIa.Xa*/
  nw_massaction_system_addReact(hm_jbc2002, 4, 5, 9, K12, K11);

  /*8*/ /* TF.VIIa + IX <-> TF.VIIa.IX -> TF.VIIa + IXa*/
  nw_massaction_system_addReact(hm_jbc2002, 4, 10, 11, K14, K13);
  /*Decay */
  nw_massaction_system_addReact(hm_jbc2002, 4, 12, 11, 0, K15);

  /*9*/ /* Xa + II -> Xa + IIa */
  /* nw_massaction_system_addReact(hm_jbc2002, 5, 13, 6, K16, 0); */
  /* nw_massaction_system_addReact(hm_jbc2002, 5, 13, 5, K16, 0); */
  nw_massaction_system_addEnzyme(hm_jbc2002, 5, 13, 6, K16);

  /*10*/ /*IIa + VIII -> IIa + VIIIa */
  /* nw_massaction_system_addReact(hm_jbc2002, 6, 14, 15, K17, 0); */
  /* nw_massaction_system_addReact(hm_jbc2002, 6, 14, 6, K17, 0); */
  nw_massaction_system_addEnzyme(hm_jbc2002, 6, 14, 15, K17);

  /*11*/ /*VIIa + IXa <-> IXa.VIIIa */
  nw_massaction_system_addReact(hm_jbc2002, 15, 12, 16, K19, K18);

  /*12*/ /*IXa.VIIIa + X <-> IXa.VIIIa.X -> IXa.VIIIa + Xa*/
  nw_massaction_system_addReact(hm_jbc2002, 16, 7, 17, K21, K20);
  /* Decay */
  nw_massaction_system_addReact(hm_jbc2002, 16, 5, 17, 0, K22);

  /*13*/ /*VIIIa1.L + VIIIa2 -> VIIIa*/
  /*Reversed...*/
  nw_massaction_system_addReact(hm_jbc2002, 18, 19, 15, K23, K24);

  /*14*/
  /* Added manually */
  /*15*/
  /* Added manually */

  /*16*/ /*IIa + V -> IIa + Va*/
  /* nw_massaction_system_addReact(hm_jbc2002, 6, 20, 21, K26, 0); */
  /* nw_massaction_system_addReact(hm_jbc2002, 6, 20, 6, K26, 0); */
  nw_massaction_system_addEnzyme(hm_jbc2002, 6, 20, 21, K26);

  /*17*/ /*Xa + Va <-> Xa.Va*/
  nw_massaction_system_addReact(hm_jbc2002, 5, 21, 22, K28, K27);

  /*18*/ /*Xa.Va + II <-> Xa.Va.II -> Xa.Va + mIIa*/
  nw_massaction_system_addReact(hm_jbc2002, 22, 13, 23, K30, K29);
  nw_massaction_system_addReact(hm_jbc2002, 22, 24, 23, 0, K31);

  /*19*/ /*mIIa + Xa.Va -> IIa + Xa.Va*/
  /* nw_massaction_system_addReact(hm_jbc2002, 24, 22, 6, K32, 0); */
  /* nw_massaction_system_addReact(hm_jbc2002, 24, 22, 22, K32, 0); */
  nw_massaction_system_addEnzyme(hm_jbc2002, 22, 24, 6, K32);

  /*20*/ /*Xa + TFPI <-> Xa.TFPI*/
  nw_massaction_system_addReact(hm_jbc2002, 5, 25, 26, K34, K33);

  /*21*/ /*TF.VIIa.Xa + TFPI <-> TF.VIIa.Xa.TFPI */
  nw_massaction_system_addReact(hm_jbc2002, 9, 25, 27, K36, K35);

  /*22*/ /*TF.VIIa + Xa.TFPI -> TF.VIIa.Xa.TFPI */
  nw_massaction_system_addReact(hm_jbc2002, 4, 26, 27, K37, 0);

  /*23*/ /*Xa + ATIII -> Xa.ATIII */
  nw_massaction_system_addReact(hm_jbc2002, 5, 28, 29, K38, 0);

  /*24*/ /*mIIa + ATIII -> mIIa.ATIII */
  nw_massaction_system_addReact(hm_jbc2002, 24, 28, 30, K39, 0);

  /*25*/ /*IXa + ATIII -> IXa.ATIII */
  nw_massaction_system_addReact(hm_jbc2002, 12, 28, 31, K40, 0);

  /*26*/ /*IIa + ATIII -> IIa.ATIII*/
  nw_massaction_system_addReact(hm_jbc2002, 6, 28, 32, K41, 0);

  /*27*/ /*TF.VIIa + TFPI -> TF.VIIa.ATIII*/
  nw_massaction_system_addReact(hm_jbc2002, 4, 28, 33, K42, 0);

  /* Add 'non-standard' reactions involving more reactants */
  /* 6 */ /*TF.VIIa.X -> TF.VIIa.Xa */
  const int R6A[]={8};
  const int R6B[]={9};
  nw_massaction_system_addGeneralReaction_S1(hm_jbc2002, R6A, 1, R6B, 1, K10, 0);

  /*14*/ /*IXa.VIIIa.X -> VIIIa1.L + VIIIa2 + X + IXa*/
  const int R14A[]={17};
  const int R14B[]={18, 19, 7, 12};
  nw_massaction_system_addGeneralReaction_S1(hm_jbc2002, R14A, 1, R14B, 4, K25, 0);

  /*15*/ /*IXa.VIIIa -> VIIIa1.L + VIIIa2 + IXa */
  const int R15A[]={16};
  const int R15B[]={18, 19, 12};
  nw_massaction_system_addGeneralReaction_S1(hm_jbc2002, R15A, 1, R15B, 3, K25, 0);

  return hm_jbc2002;
}

#ifdef TEST
int main(int argc, char **argv) {

  nw_massaction_system *  hm_jbc2002=init_hm_jbc2002();
  /* nw_massaction_system_alloc(NCOMP, NREACT); */

  gsl_odeiv2_system gslModel = {nw_massaction_system_ODEs/* hm2002 */,
				/* jac_hm2002jac */nw_massaction_system_jac, NCOMP, hm_jbc2002};
  gsl_odeiv2_driver * d;
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

  int i, j;
  double t=0, t1=100, tend=/* 50000 */700/* 22279 *//* 10000 */;
  int nBlood=34;
  gsl_matrix *Smatrix;
  int withTFPI=0, withATIII=0;

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
  d = gsl_odeiv2_driver_alloc_y_new (&gslModel, gsl_odeiv2_step_rkf45/* gsl_odeiv2_step_msbdf *//* gsl_odeiv2_step_bsimp *//* gsl_odeiv2_step_rk8pd */,
				     1e-6, 1e-6, 0.0);

  /* Title */
  fprintf(OUT, "t\t");
  for(i=0; i<nBlood-1; i++) fprintf(OUT, "%s\t", nw_massaction_system_getName(hm_jbc2002, i));
  fprintf(OUT, "%s\n", nw_massaction_system_getName(hm_jbc2002, i));

  fprintf(OUT, "%g\t", 0.0);
  for(i=0; i<nBlood-1; i++) fprintf(OUT, "%g\t", y[i]);
  fprintf(OUT, "%g\n", y[i]);

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
    }

  gsl_odeiv2_driver_free (d);
  

  return 0;  
}
#endif
