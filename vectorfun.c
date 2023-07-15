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

const float PI=3.141593;

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

double arccos(double x) {
  return acos(x)*180.0/PI;
}

double get_angle(double *v, double* w) {
  return acos(dot(v,w)/(norm(v)*norm(w)))*180.0/PI;
}

double argument(double x, double y) {
  return atan(y/x)*180.0/PI;
}

double* orthonormal_projection(double *v, double *w) {
  double* result = malloc(3 * sizeof(double));
  double v_dot_w = dot(v,w);
  double w_dot_w = dot(w,w);
  result[0] = v[0] - w[0] * v_dot_w / w_dot_w;
  result[1] = v[1] - w[1] * v_dot_w / w_dot_w;
  result[2] = v[2] - w[2] * v_dot_w / w_dot_w;
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
