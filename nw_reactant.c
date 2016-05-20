/*
 *  nw_reactant.c
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nw_reactant.h"

#define OUT stdout
#define ERR stderr

/* Allocate memory for reactant */
nw_reactant *nw_reactant_alloc(char * name, int ri, int s) {
  nw_reactant *r=malloc(sizeof(*r));

  nw_reactant_setR(r, ri);
  nw_reactant_setStoichiometry(r, s);
  nw_reactant_setName(r, name);
  return r;
}

/* Sets name */
void nw_reactant_setName(nw_reactant *r, char *name) {
  r->name=name;
}

/* Returns name of reactant */
char *nw_reactant_getName(const nw_reactant *r) {
  return r->name;
}

/* Return index of reactant */
int nw_reactant_getR(const nw_reactant *r) {
  return r->r;
}

/* Set index of reactant */
void nw_reactant_setR(nw_reactant *r, int ri) {
  r->r=ri;
}

/* Set index of reactant with stoichiometry 1.*/
void nw_reactant_setSingleR(nw_reactant *r, int ri) {
  nw_reactant_setR(r, ri);
  nw_reactant_setStoichiometry(r, 1);
}


/* Return stoichiometry... as int */
int nw_reactant_getStoichiometry(const nw_reactant *r) {
  return (int)r->s;
}

/* Set stoichiometry... as int */
void nw_reactant_setStoichiometry(nw_reactant *r, int s) {
  r->s=s;
}

/* String representation of reactant */
char *nw_reactant_toString(const nw_reactant *r) {
  char *name=nw_reactant_getName(r);
  char *out=malloc(3+strlen(name) + 1);
  int s=nw_reactant_getStoichiometry(r);

  sprintf(out, "%2i %s", s, name);
  return out;
}
