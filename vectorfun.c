/*
 * very simple library for some basic vector operations. 
 * only for 3-dimensional vectors
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#include "vectorfun.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

double* add(double* v, double* w) {
  double* result = malloc(3 * sizeof(double));
  result[0] = v[0] + w[0];
  result[1] = v[1] + w[1];
  result[2] = v[2] + w[2];
  return result;
}

double* subtract(double* v, double* w) {
  double* result = malloc(3 * sizeof(double));
  result[0] = v[0] - w[0];
  result[1] = v[1] - w[1];
  result[2] = v[2] - w[2];
  return result;
}

double dot(double* v, double* w) {
  assert(v && w);
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

double norm(double *v) {
  assert(v);
  return sqrt(dot(v,v));
}

double* normalize(double *v) {
  assert(v);
  double* result = malloc(3 * sizeof(double));
  result[0] = v[0]/norm(v);
  result[1] = v[1]/norm(v);
  result[2] = v[2]/norm(v);
  return result;
}

double get_angle(double *v, double* w) {
  printf("dot: %f\n",dot(v,w));
  return acos(dot(v,w));
}

double* orthonormal_projection(double *v, double *w) {
  double* result = malloc(3 * sizeof(double));
  double dot_product = dot(v,w);
  double norm_w = norm(w);
  result[0] = v[0] - w[0] * dot_product / dot(w,w);
  result[1] = v[1] - w[1] * dot_product / dot(w,w);
  result[2] = v[2] - w[2] * dot_product / dot(w,w);
  return normalize(result);
}

double* cross(double* v, double* w) {
  assert(v && w);
  double* result = malloc(3 * sizeof(double));
  result[0] = v[1]*w[2] - v[2]*w[1];
  result[1] = v[2]*w[0] - v[0]*w[2];
  result[2] = v[0]*w[1] - v[1]*w[0];
  return result;
}
