/*
 * very simple library for some basic vector operations. 
 * only for 3-dimensional vectors
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#ifndef _VECTORFUN
#define _VECTORFUN

#include <math.h>
#include <stdlib.h>
#include <assert.h>

double* add(double* v, double* w);
double* subtract(double* v, double* w);
double dot(double* v, double* w);
double norm(double *v);
double* normalize(double *v);
double arccos(double x);
double get_angle(double *v, double* w);
double argument(double x, double y);
double* orthonormal_projection(double *v, double *w);
double* cross(double* v, double* w);
double* matvec_product(double** A, double* x);
double** matmat_product(double** A, double** B);
double** inverse(double** A);

#endif
