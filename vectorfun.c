/*
 * very simple library for some basic vector operations. 
 * only for 3-dimensional vectors
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#include <stdlib.h>
#include <assert.h>

double dot_product(double* v, double* w) {
  assert(v && w);
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

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

double* cross_product(double* v, double* w) {
  assert(v && w);
  double* result = malloc(3 * sizeof(double));
  result[0] = v[1]*w[2] - v[2]*w[1];
  result[1] = v[2]*w[0] - v[0]*w[2];
  result[2] = v[0]*w[1] - v[1]*w[0];
  return result;
}
