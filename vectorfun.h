/*
 * very simple library for some basic vector operations. 
 * only for 3-dimensional vectors
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 */

#ifndef _VECTORFUN
#define _VECTORFUN

#include <math.h>
#include <stdlib.h>
#include <assert.h>

extern const double pi;

// vector methods

double* create_vector();
void free_vector(double* v);
void set_vector(double* v, double* w);
double* copy_vector(double* v);
double* add(double* v, double* w);
double* subtract(double* v, double* w);
double* scale_vector(double *v, double scalar);
double dot(double* v, double* w);
double norm(double *v);
double* normalize(double *v);
double* orthonormal_projection(double *v, double *w);
double* cross(double* v, double* w);
void print_vector(double* v);

// matrix methods 

double** create_matrix(double val);
void free_matrix(double** A);
void set_matrix(double** A, double** B);
double** copy_matrix(double** A);
void set_matrix_rows(double** A, double* x, double* y, double* z);
double* matvec_product(double** A, double* x);
double** matmat_product(double** A, double** B);
double** inverse(double** A);
double** transpose(double** A);
void print_matrix(double** A);

// angle methods 

double arccos(double x);
double get_angle(double *v, double* w);
double argument(double x, double y);
double* compute_euler_angles(double** R);
double* compute_euler_vector(double** R);
double* compute_euler_vector2(double** R);

#endif
