#ifndef NW_REACTANT
#define NW_REACTANT
/*
 *  nw_reactant.h
 *
 *
 * Type for reactant in a chemical reaction
 * 
 *
 * Ivo Siekmann, 11.04.2014
 *
 * 
 * Compile with:
 *
 * gcc -Wall -pedantic -o filename filename.c
 *
 */

typedef struct {
  /* Name of reactant */
  char * name;
  /* Index of reactant */
  int r;
  /* Stoichometry: Double for including "fractional" stoichiometries */
  double s;
} nw_reactant;

/* Allocate memory for reactant */
nw_reactant *nw_reactant_alloc(char *name, int ri, int s);

/* Sets name */
void nw_reactant_setName(nw_reactant *r, char *name);

/* Returns name of reactant */
char *nw_reactant_getName(const nw_reactant *r);

/* Return index of reactant */
int nw_reactant_getR(const nw_reactant *r);

/* Set index of reactant */
void nw_reactant_setR(nw_reactant *r, int ri);

/* Set index of reactant with stoichiometry 1.*/
void nw_reactant_setSingleR(nw_reactant *r, int ri);

/* Return stoichiometry... as int */
int nw_reactant_getStoichiometry(const nw_reactant *r);

/* Set stoichiometry... as int */
void nw_reactant_setStoichiometry(nw_reactant *r, int s);

/* String representation of reactant */
char *nw_reactant_toString(const nw_reactant *r);

#endif
