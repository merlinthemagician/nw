#ifndef JONESMANN
#define JONESMANN
/*
 *  JonesMann_JBC1994.h
 *
 *
 * Jones & Mann coagulation model
 *
 * Model parameters
 *
 * Ivo Siekmann, 07/03/2014
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

#define NCOMP 30
#define NREACT 40

#define SCL (1e-9)

#define K1 (2e7*SCL)
#define K2 (2e7*SCL)
#define K3 (1e7*SCL)
#define K4 (2e7*SCL)
#define K5 (1e7*SCL)
#define K6 (1e8*SCL)
#define K7 (1e7*SCL)
#define K8 (4e8*SCL)

#define K9 (0.005)
#define K10 (0.4)
#define K11 (0.3)
#define K12 (1.15)
#define K13 (8.2)
#define K14 (32)

#define K15 (1e5*SCL)
#define K16 (24)
#define K17 (44)
#define K18 (0.001)
#define K19 (70)
#define K20 (0.02)

/* Plasma concentrations [nM] */
#define II0 (1400.0)
#define V0 (20.0)
#define VIII0 (0.7)
#define VII0 (10.0)
#define IX0 (90.0)
#define X0 (170.0)
#define XI0 (31.2)

const char *compNames[] = { "TF.VIIa", 
			    "IX",
			    "X",
			    "V",
			    "VIII",
			    "II",
			    "VIIIa.IXa",
			    "Va.Xa",
			    "IIa",
			    "Va.Xa.II",
			    "mIIa",
			    "TF.VIIa.IX",
			    "TF.VIIa.X",
			    "VIIIa.IXa.X",
			    "IXa",
			    "Xa",
			    "Va",
			    "VIIIa",
			    "I"
};

#endif
