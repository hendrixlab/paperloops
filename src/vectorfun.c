/*
 * very simple library for some basic vector operations. 
 * only for 3-dimensional vectors
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 */

#include "vectorfun.h"
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

const double pi = 3.1415926535898;

/*******************
 * General Methods *
 *******************/

double max(double x, double y) {
  if(x >= y) {
    return x;
  } else {
    return y;
  }
}

/******************
 * Vector Methods *
 *******************/

double* create_vector() {
  double* vector = malloc(3 * sizeof(double));
  vector[0] = 0.0;
  vector[1] = 0.0;
  vector[2] = 0.0;
  return vector;
}

void free_vector(double* v) {
  assert(v);
  free(v);
}

void set_vector(double* v, double* w) {
  assert(v);
  assert(w);
  v[0] = w[0];
  v[1] = w[1];
  v[2] = w[2];
}

void set_vector_elements(double* v, double x, double y, double z) {
  assert(v);
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

double* copy_vector(double* v) {
  double* copy = create_vector();
  copy[0] = v[0];
  copy[1] = v[1];
  copy[2] = v[2];
  return copy;
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

void print_vector(double* v ) {
  printf("[");
  for(int i=0;i<3;i++) { 
    printf(" %.3f ", v[i]);
  }
  printf("]\n");
}

/******************
 * Matrix Methods *
 ******************/

double** create_matrix(double val) {
  // when val is 0.0 it is a zero matrix, when 1.0 it's an identity matrix
  double** matrix = malloc(3 * sizeof(double*));
  matrix[0] = malloc(3*sizeof(double));
  matrix[1] = malloc(3*sizeof(double));
  matrix[2] = malloc(3*sizeof(double));
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      if(i == j) {
	matrix[i][j] = val;
      } else {
	matrix[i][j] = 0.0;
      }
    }
  }
  return matrix;
}

void free_matrix(double** A) {
  assert(A);
  assert(A[0]);
  assert(A[1]);
  assert(A[2]);
  free(A[0]);
  free(A[1]);
  free(A[2]);
  free(A);
}

void set_matrix(double** A, double** B) {
  assert(A);
  assert(B);
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      A[i][j] = B[i][j];
    }
  }
}

double** copy_matrix(double** A) {
  double** copy = create_matrix(0.0);
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      copy[i][j] = A[i][j];
    }
  }
  return copy;
}

 void set_matrix_rows(double** A, double* x, double* y, double* z) {
  for(int j=0; j<3; j++) {
    A[0][j] = x[j];
    A[1][j] = y[j];
    A[2][j] = z[j];
  }
}

double* matvec_product(double** A, double* x) {
  double* result = create_vector();  
  result[0] = dot(A[0],x);
  result[1] = dot(A[1],x);
  result[2] = dot(A[2],x);
  return result;
}

double** matmat_product(double** A, double** B) {
  double** result = create_matrix(0.0);
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
  double** inv = create_matrix(1.0);
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
  double** A_Tr = create_matrix(1.0);
  
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++) {
      A_Tr[i][j] = A[j][i];
    }  
  }
  return A_Tr;
}

double trace(double** A) {
  // returns the trace of a matrix (sum of diagonal terms)
  double Tr = 0.0;
  for(int i=0; i<3; i++) {
    Tr += A[i][i];
  }
  return Tr;
}

void print_matrix(double** A) {
  for(int i=0;i<3;i++) {
    printf("[");
    for(int j=0;j<3;j++) {
      printf(" %.3f ",A[i][j]);
    }
    printf("]\n");
  }
}

/*****************
 * Angle Methods *
 *****************/

double arccos(double x) {
  return acos(x)*180.0/pi;
}

double get_angle(double *v, double* w) {
  return acos(dot(v,w)/(norm(v)*norm(w)))*180.0/pi;
}

double argument(double x, double y) {
  return atan(y/x)*180.0/pi;
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

/*
This function is no longer used: using acos() left the angle ambiguous by a factor of pi
Replaced below using intermediate quaternion calculation that solves this. 

double* compute_euler_vector(double** R) {
  // compute the Euler vector from rotation matrix, i.e. to "Euler axis/angle"
  // defined as the product of the axis and angle from axis-angle representation
  // clearest formula: https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions
  // FLAW: assumes that theta is not an exact mulitple of pi, maybe unstable around there
  double theta = acos((trace(R)-1.0)/2.0);
  //printf("Tr(R)=%.3f, theta=%.3f\n",trace(R),theta);
  double* axis = create_vector(); // allocate 0-vector
  axis[0] = (R[2][1] - R[1][2])/(2.0*sin(theta));
  axis[1] = (R[0][2] - R[2][0])/(2.0*sin(theta));
  axis[2] = (R[1][0] - R[0][1])/(2.0*sin(theta));  
  double* EV = scale_vector(axis,theta);
  //print_vector(EV);
  return EV;
}
*/

double* compute_euler_vector(double** R) {
  // compute the Euler vector from rotation matrix, i.e. to "Euler axis/angle"
  // converts to quaternion first, based on: https://intra.engr.ucr.edu/~farrell/AidedNavigation/D_App_Quaternions/Rot2Quat.pdf
  double* q = malloc(4*sizeof(double));
  q[0] = 0.5*sqrt(1.0+R[0][0]-R[1][1]-R[2][2]); // q_i
  q[1] = 0.5*sqrt(1.0-R[0][0]+R[1][1]-R[2][2]); // q_j
  q[2] = 0.5*sqrt(1.0-R[0][0]-R[1][1]+R[2][2]); // q_k
  q[3] = 0.5*sqrt(1.0+R[0][0]+R[1][1]+R[2][2]); // q_r the angle/real component
  // get the largest term to maximize numerical precision
  int max_index = 0;
  for(int index=1; index<4; index++) {
    if(q[max_index] < q[index]) {
      max_index = index;
    }
  }
  // update other terms using the max denominator
  switch(max_index) {
  case 0:
    q[1] = (R[0][1] + R[1][0])/(4.0*q[0]);
    q[2] = (R[0][2] + R[2][0])/(4.0*q[0]);
    q[3] = (R[1][2] - R[2][1])/(4.0*q[0]);
  case 1:
    q[0] = (R[0][1] + R[1][0])/(4.0*q[1]);
    q[2] = (R[1][2] + R[2][1])/(4.0*q[1]);
    q[3] = (R[2][0] - R[0][2])/(4.0*q[1]);
  case 2:
    q[0] = (R[0][2] + R[2][0])/(4.0*q[2]);
    q[1] = (R[1][2] + R[2][1])/(4.0*q[2]);
    q[3] = (R[0][1] - R[1][0])/(4.0*q[2]);
  case 3:
    q[0] = (R[1][2] - R[2][1])/(4.0*q[3]);
    q[1] = (R[2][0] - R[0][2])/(4.0*q[3]);
    q[2] = (R[0][1] - R[1][0])/(4.0*q[3]);
  }
  // compute axis-angle representation from quaternions
  // based on: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
  double norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
  double theta = 2.0*atan2(norm,q[3]);
  double* axis = create_vector();
  if(sin(theta/2.0) == 0.0) {
    axis[0] = 1.0; // the other terms are 0.0 by default
  } else {
    axis[0] = q[0]/sin(theta/2.0);
    axis[1] = q[1]/sin(theta/2.0);
    axis[2] = q[2]/sin(theta/2.0);
  }
  //  printf("%f\t%f\t%f\t%f\n",axis[0],axis[1],axis[2],theta);
  double* EV = scale_vector(axis,theta);
  // free memory
  free(q);
  free(axis);
  return EV;
}

