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

const double pi = 3.1415926535898;

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

double* scale_vector(double *v, double scale) {
  // multiply a vector by a scalar
  double* result = create_vector();
  result[0] = scale*v[0];
  result[1] = scale*v[1];
  result[2] = scale*v[2];
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
  return acos(x)*180.0/pi;
}

double get_angle(double *v, double* w) {
  return acos(dot(v,w)/(norm(v)*norm(w)))*180.0/pi;
}

double argument(double x, double y) {
  return atan(y/x)*180.0/pi;
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

/* Matrices */

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
  
  for(int i=0; i<3; i++) {
    det = det + (A[0][i]*(A[1][(i+1)%3]*A[2][(i+2)%3] - A[1][(i+2)%3]*A[2][(i+1)%3]));
  }  
  for(int i=0; i<3; i++){
    for(int j=0;j<3;j++) {
      inv[i][j] = (A[(j+1)%3][(i+1)%3]*A[(j+2)%3][(i+2)%3] - A[(j+2)%3][(i+1)%3]*A[(j+1)%3][(i+2)%3]) / det;
    }  
  }
  return inv;
}

double** transpose(double** A) {
  double** A_Tr = create_matrix();
  
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++) {
      A_Tr[i][j] = A[j][i];
    }  
  }
  return A_Tr;
}

double trace(double** A) {
  // returns the trace of a matrix (sum of diagonal terms)
  double Tr = 0;
  for(int i=0; i<3; i++) {
    Tr += A[i][i];
  }
  return Tr;
}

double* compute_euler_angles(double** R) {
  // compute Euler angles from rotation matrix
  // based on math from: https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf
  double chi = atan2(-R[2][0],sqrt(R[0][0]*R[0][0]+R[1][0]*R[1][0]));
  double phi = atan2(R[1][0],R[0][0]);
  double psi = atan2(R[2][1],R[2][2]);
  if(chi == -90.0) {
    phi = atan2(-R[1][2],-R[0][2]);
    psi = 0.0;
  } else if(chi == 90.0) {
    phi = atan2(R[1][2],R[0][2]);
    psi = 0.0;
  }
  double* EA = create_vector();
  EA[0] = chi;
  EA[1] = phi;
  EA[2] = psi;
  return EA;
}

double* compute_euler_vector(double** R) {
  // compute the Euler vector from rotation matrix $
  // defined as the product of the axis and angle from axis-angle representation
  // clearest formula: https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions
  // FLAW: assumes that theta is not an exact mulitple of pi, maybe unstable around there
  double theta = acos((trace(R)-1.0)/2.0);
  double* omega = create_vector(); // allocate 0-vector
  omega[0] = (R[2][1] - R[1][2])/(2.0*sin(theta));
  omega[1] = (R[0][2] - R[2][0])/(2.0*sin(theta));
  omega[2] = (R[1][0] - R[0][1])/(2.0*sin(theta));
  double* EV = scale_vector(omega,theta);
  return EV;
}

