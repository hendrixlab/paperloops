/*
 * very simple library for some basic vector operations. 
 * only for 3-dimensional vectors
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#include <stdlib.h>
#include <assert.h>

double dot_product(double* v, double* w);
double* add(double* v, double* w);
double* subtract(double* v, double* w);
double* cross_product(double* v, double* w);
