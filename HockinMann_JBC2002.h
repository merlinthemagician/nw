#ifndef HOCKINMANN
#define HOCKINMANN
/*
 *  HockinMann_JBC1994.h
 *
 *
 * Hockin et al. & Mann coagulation model
 *
 * Model parameters
 *
 * Ivo Siekmann, 17/03/2014
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#define NCOMP 34
#define NREACT 40

#define SCL (1e-9)

#define K1 (3.1e-3)

#define K2 (3.2e6*SCL)

#define K3 (3.1e-3)

#define K4 (2.3e7*SCL)

/*Corrected from paper */
#define K5 (4.4e5*SCL)
#define K6 (1.3e7*SCL)
#define K7 (2.3e4*SCL)

#define K8 (1.05)

#define K9 (2.5e7*SCL)

#define K10 (6)
#define K11 (19)

#define K12 (2.2e7*SCL)

#define K13 (2.4)

#define K14 (1e7*SCL)

#define K15 (1.8)

#define K16 (7.5e3*SCL)
#define K17 (2e7*SCL)

#define K18 (5e-3)

#define K19 (1e7*SCL)

#define K20 (1e-3)

#define K21 (1e8*SCL)

#define K22 (8.2)

#define K23 (2.2e4*SCL)

#define K24 (6e-3)
#define K25 (1e-3)

#define K26 (2e7*SCL)

#define K27 (0.2)

#define K28 (4e8*SCL)

#define K29 (103)

#define K30 (1e8*SCL)

#define K31 (63.5)

#define K32 (1.5e7*SCL)

#define K33 (3.6e-4)

#define K34 (9e5*SCL)

#define K35 (1.1e-4)

#define K36 (3.2e8*SCL)
#define K37 (5e7*SCL)
#define K38 (1.5e3*SCL)
#define K39 (7.1e3*SCL)
#define K40 (4.9e2*SCL)
#define K41 (7.1e3*SCL)
#define K42 (2.3e2*SCL)

/* Plasma concentrations [nM] */
#define VII0 (10.0)
#define VIIa0 (0.1)

#define X0 (160.0)

#define IX0 (90.0)

#define II0 (1400.0)

#define VIII0 (0.7)

#define V0 (20.0)

#define TFPI0 (2.5)

#define ATIII0 (3400)

static const char *compNames[] = { "TF", /*0*/
			    "VII", /*1*/ 
			    "TF.VII", /*2*/
			    "VIIa", /*3*/
			    "TF.VIIa", /*4*/
			    "Xa", /*5*/
			    "IIa", /*6*/
			    "X", /*7*/
			    "TF.VIIa.X", /*8*/
			    "TF.VIIa.Xa", /*9*/
			    "IX",/*10*/
			    "TF.VIIa.IX", /*11*/
			    "IXa", /*12*/
			    "II", /*13*/
			    "VIII", /*14*/
			    "VIIIa", /*15*/
			    "IXa.VIIIa", /*16*/
			    "IXa.VIIIa.X", /*17*/
			    "VIIIa1.L", /*18*/
			    "VIIIa2", /*19*/ 
			    "V", /*20*/
			    "Va", /*21*/
			    "Xa.Va", /*22*/ 
			    "Xa.Va.II", /*23*/
			    "mIIa", /*24*/
			    "TFPI", /*25*/
			    "Xa.TFPI", /*26*/
			    "TF.VIIa.Xa.TFPI", /*27*/
			    "ATIII", /*28*/
			    "Xa.ATIII", /*29*/
			    "mIIa.ATIII", /*30*/
			    "IXa.ATIII", /*31*/
			    "IIa.ATIII", /*32*/
			    "TF.VIIa.ATIII" /*33*/
};

/* conserved quantities */
static int C1[]={0, 2, 4, 8, 9, 11, 27, 33};
static int C2[]={1, 2, 3, 4, 8, 9, 11, 27, 33};
static int C3[]={5, 7, 8, 9, 17, 22, 23, 26, 27, 29};
static int C4[]={6, 13, 23, 24, 30, 32};
static int C5[]={10, 11, 12, 16, 17, 31};
static int C6[]={14, 15, 16, 17, 19};
/* Equality. NOT sum!*/
static int C7[]={18, 19};
static int C8[]={20, 21, 22, 23};
static int C9[]={25, 26, 27};
static int C10[]={28, 29, 30, 31, 32, 33};

static double kerHM[] = {
  1,0,1,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,
  0,1,1,1,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,
  0,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,0,1,0,0,0,0,
  0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,0,
  0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1
};

/* Differential equations of the Hockin-Mann model: par is expected to be an nw_massaction_system */
int hm2002(double t, const double y[], double dydt[], void *par);

/* Jacobian of the Hockin-Mann model: par is expected to be an nw_massaction_system */
int jac_hm2002jac (double t, const double y[], double *dfdy, 
		   double dfdt[], void *par);

/* Initialise Hockin/Mann system */
nw_massaction_system* init_hm_jbc2002();
#endif
