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

const float PI=3.14159265359;

double* create_vector() {
  double* vector = malloc(3 * sizeof(double));
  vector[0] = 0.0;
  vector[1] = 0.0;
  vector[2] = 0.0;
  return vector;
}

double** create_matrix() {
  double** matrix = malloc(3 * sizeof(double*));
  matrix[0] = malloc(3*sizeof(double));
  matrix[1] = malloc(3*sizeof(double));
  matrix[2] = malloc(3*sizeof(double));
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      if(i == j) {
	matrix[i][j] = 1.0;
      } else {
	matrix[i][j] = 0.0;
      }
    }
  }
  return matrix;
}

double* add(double* v, double* w) {
  double* result = create_vector();
  result[0] = v[0] + w[0];
  result[1] = v[1] + w[1];
  result[2] = v[2] + w[2];
  return result;
}

double* subtract(double* v, double* w) {
  double* result = create_vector();  
  result[0] = v[0] - w[0];
  result[1] = v[1] - w[1];
  result[2] = v[2] - w[2];
  return result;
}

double* negate(double *v) {
  double* result = create_vector();
  result[0] = -v[0];
  result[1] = -v[1];
  result[2] = -v[2];
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
  double* result = create_vector();
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
  double* result = create_vector();  
  double v_dot_w = dot(v,w);
  double w_dot_w = dot(w,w);
  result[0] = v[0] - w[0] * v_dot_w / w_dot_w;
  result[1] = v[1] - w[1] * v_dot_w / w_dot_w;
  result[2] = v[2] - w[2] * v_dot_w / w_dot_w;
  return normalize(result);
}

double* cross(double* v, double* w) {
  assert(v && w);
  double* result = create_vector();  
  result[0] = v[1]*w[2] - v[2]*w[1];
  result[1] = v[2]*w[0] - v[0]*w[2];
  result[2] = v[0]*w[1] - v[1]*w[0];
  return result;
}

double* matvec_product(double** A, double* x) {
  double* result = create_vector();  
  result[0] = dot(A[0],x);
  result[1] = dot(A[1],x);
  result[2] = dot(A[2],x);
  return result;
}

double** matmat_product(double** A, double** B) {
  double** result = create_matrix();
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      for(int k=0; k<3; k++) {
	result[i][j] += A[i][k]*B[k][j];
      }
    }
  }
  return result;
}

double** inverse(double** A) {
  double** inv = create_matrix();
  double det = 0.0;
  
  for(int i=0;i<3;i++) {
    det = det + (A[0][i]*(A[1][(i+1)%3]*A[2][(i+2)%3] - A[1][(i+2)%3]*A[2][(i+1)%3]));
  }  
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) {
      inv[i][j] = (A[(j+1)%3][(i+1)%3]*A[(j+2)%3][(i+2)%3] - A[(j+2)%3][(i+1)%3]*A[(j+1)%3][(i+2)%3]) / det;
    }  
  }
  return inv;
}

double** transpose(double** A) {
  double** A_Tr = create_matrix();
  
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) {
      A_Tr[i][j] = A[j][i];
    }  
  }
  return A_Tr;
}

