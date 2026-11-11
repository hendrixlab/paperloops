 /*
 * This library contains the definitions of structures and functions for transforms between 
 * adjacent nucleotides
 *
 * (C) David Hendrix david.hendrix@oregonstate.edu
 * In collaboration with Rhiju Das, Stanford
 * translated/derived from toyfold3 with permission: https://github.com/rhiju/toyfold3 (C) Rhiju Das
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "transform.h"
#include "rna_data.h"
#include "motifs.h"
#include "vectorfun.h"

//const double pi = 3.1415926535898;

/**********************
 * Methods for frames *
 **********************/

struct frame* create_frame(double val) {
  struct frame* f = malloc(sizeof(struct frame));
  f->ctr = create_vector(); // vector of zeros, placeholer
  f->M = create_matrix(val); // val along the diagonal, placeholder
  f->EV = create_vector(); // Euler vector, initially vector of zeros
  return f;
}

void free_frame(struct frame* f) {
  assert(f);
  free_vector(f->ctr);
  free_matrix(f->M);
  free_vector(f->EV);
  free(f);
}

struct frame* compute_frame(struct nuc* n) {
  // create a frame
  struct frame* f = create_frame(1.0);
  // compute frame translation and rotation
  struct atom* a = get_backbone_atom_by_name(n,"C5'"); 
  struct atom* b = get_backbone_atom_by_name(n,"C4'"); 
  struct atom* c = get_backbone_atom_by_name(n,"C3'");
  // compute x,y, and z axes from a,b,c
  double* b_minus_a = subtract(atom_coord(b),atom_coord(a));
  double* z = normalize(b_minus_a);
  double* c_minus_b = subtract(atom_coord(c),atom_coord(b));
  double* cross_prod1 = cross(z,c_minus_b);
  double* y = normalize(cross_prod1);
  double* cross_prod2 = cross(y,z);
  double* x = normalize(cross_prod2);
  // build frame orientation matrix from x, y and z
  double** M = create_matrix(1.0);
  set_matrix_rows(M,x,y,z);
  double** M_T = transpose(M);
  // define the frame values:
  set_vector(f->ctr,atom_coord(a));
  set_matrix(f->M,M_T);
  set_frame_num(f,nuc_num(n));
  // some memory cleanup:
  free_vector(b_minus_a);
  free_vector(c_minus_b);
  free_vector(cross_prod1);
  free_vector(cross_prod2);
  free_vector(x);
  free_vector(y);
  free_vector(z);
  free_matrix(M);
  free_matrix(M_T);
  return f;
}

struct frame* reverse_frame(struct frame* f) {
  struct frame* f_rev = create_frame(1.0);
  double** M_inv = inverse(f->M);  // transpose should work and seems faster
  double* ctr_neg = scale_vector(f->ctr,-1.0);
  double* ctr_rev = matvec_product(M_inv,ctr_neg);
  set_vector(f_rev->ctr,ctr_rev);
  set_matrix(f_rev->M,M_inv);
  free_vector(ctr_neg);
  free_vector(ctr_rev);
  free_matrix(M_inv);
  return f_rev;
}

struct frame* copy_frame(struct frame* f) {
  struct frame* f_copy = create_frame(1.0);
  set_vector(f_copy->ctr,f->ctr);
  set_matrix(f_copy->M,f->M);
  set_vector(f_copy->EV,f->EV);
  return f_copy;
}

struct frame* average_frames(struct frame* f1, struct frame* f2) {
  struct frame* f_avg = create_frame(1.0);
  double* ctr1 = get_frame_center(f1);
  double* ctr2 = get_frame_center(f2);
  for(int i=0; i<3; i++) {
    f_avg->ctr[i] = (ctr1[i]+ctr2[i])/2.0;
  }
  return f_avg;
}

double get_frame_distance_diff(struct frame* f1, struct frame* f2) {
  double* diff = subtract(get_frame_center(f1),get_frame_center(f2));
  double distance = norm(diff);
  free_vector(diff);
  return distance;
}

double get_frame_angle_diff(struct frame* f1, struct frame* f2) {
  double** M1_inv = inverse(get_frame_orientation(f1));
  double** R_diff = matmat_product(M1_inv,get_frame_orientation(f2));
  double* EV = compute_euler_vector(R_diff);
  double angle = norm(EV);
  free_matrix(M1_inv);
  free_matrix(R_diff);
  free_vector(EV);
  return angle;
}

double* get_frame_center(struct frame* f) {
  return f->ctr;
}

void set_frame_center(struct frame* f, double* ctr) {
  set_vector(f->ctr,ctr);
}

double** get_frame_orientation(struct frame* f) {
  return f->M;
}

void set_frame_orientation(struct frame* f, double** M) {
  set_matrix(f->M,M);
}

void add_frame_euler_vector(struct frame* f) {
  double* EV = compute_euler_vector(get_frame_orientation(f));
  set_frame_euler_vector(f,EV);
  free_vector(EV); // set_vector assigns values element-by-element so original vector needs to be freed
}

void set_frame_euler_vector(struct frame* f, double* EV) {
  //printf("test: %.3f %.3f %.3f\n",EV[0],EV[1],EV[2]);
  set_vector(f->EV,EV);
}

double* get_frame_euler_vector(struct frame* f) {
  return f->EV;
}

int get_frame_num(struct frame* f) {
  return f->num;
}

void set_frame_num(struct frame* f, int num) {
  f->num = num;
}

// T6 Methods

double get_T6_value(struct frame* f, int i) {
  if(0 <= i && i < 3) {
    return f->ctr[i];
  } else if(3 <= i && i < 6) { // hard-coded b/c should never be other than 6D
    return f->EV[i-3];
  } else {
    printf("Unexpected value in get_T6_value: %d",i);
    exit(0);
  }
}

double* get_T6_std_dev(struct trace* pts) {
  int d = 6; // assuming 6D coordinates
  double* sum = malloc(sizeof(double)*d);
  double* sum_squares = malloc(sizeof(double)*d);
  double* std_dev = malloc(sizeof(double)*d);  
  int n = trace_length(pts);
  
  for(int j=0; j<d; j++) {
    sum[j] = 0.0;
    sum_squares[j] = 0.0;
  }
  for(int i=0; i<n; i++) {
    struct frame* f = get_frame(pts,i);
    for(int j=0; j<d; j++) {
      double val = get_T6_value(f,j);
      sum[j] += val;
      sum_squares[j] += val*val;	
    }
  } 
  for(int j=0; j<d; j++) {
    std_dev[j] = sqrt((sum_squares[j] - sum[j]*sum[j]/n)/(n-1));
  }
  free(sum);
  free(sum_squares);
  
  return std_dev;
}

void print_frame_T6(struct frame* f) {
  int d=6;
  printf("[");
  for(int j=0; j<d; j++) {
    double x = get_T6_value(f,j);
    printf(" %f ",x);
  }
  printf("]\n");
}


/**************************************
 * Tools for traces (lists of frames) *
 **************************************/

struct trace* create_trace(int length) {
  struct trace* trace = malloc(sizeof(struct trace));
  assert(trace);

  trace->data = malloc(length*sizeof(struct frame*));
  assert(trace->data);
  trace->len = length;

  return trace;
}

void free_trace(struct trace* trace) {
  assert(trace);  
  for(int i=0; i<trace_length(trace); i++) {
    struct frame* f = get_frame(trace,i);
    free_frame(f);
  }
  free(trace->data);
  free(trace);
}

int trace_length(struct trace* trace) {
  assert(trace);
  return trace->len;
}
      
struct frame* get_frame(struct trace* trace, int idx) {
  assert(trace);
  if(!((0 <= idx && idx < trace->len) || idx == -1)) {
    printf("get_frame: attempted idx=%d vs trace_length=%d\n",idx,trace->len);
    exit(1);
  }

  if (idx == -1) { //  idx = -1 grabs from the end of the array.
    idx = trace->len - 1;
  }

  return trace->data[idx];
}

void set_frame(struct trace* trace, int idx, struct frame* f) {
  assert(trace);

  if((0 <= idx)&&(idx < trace_length(trace))) {
    trace->data[idx] = f;
  } else {
    printf("index out of bounds. Attempted %d for trace of size %d\n",idx,trace_length(trace));
    exit(1);
  }
}

double* get_trace_min_bounds(struct trace* trace, double* input_min) {
  // min vector passed in and updated to allow custom initial value for min
  assert(trace);
  // initialization, optional input_min for running min updates
  double* min = create_vector();
  if(input_min) {
    min = input_min;
  } else {
    min[0] = min[1] = min[2] = 1000000000.0;    
  }
  // find actual min from the trace
  for(int i=0; i<trace_length(trace); i++) {
    struct frame* f = get_frame(trace,i);
    double* x = get_frame_center(f);
    //printf("x=(%.2f,%.2f,%.2f)\n",x[0],x[1],x[2]);    
    for(int d=0; d<3; d++) {
      if(x[d] < min[d]) {
	min[d] = x[d];
      }
    }
  }
  return min;
}

double* get_trace_max_bounds(struct trace* trace, double* input_max) {
  assert(trace);
  // set up initial values
  double* max = create_vector();
  if(input_max) {
    max = input_max;
  } else {
    max[0] = max[1] = max[2] = -1000000000.0;
  }
  // find max from trace
  for(int i=0; i<trace_length(trace); i++) {
    struct frame* f = get_frame(trace,i);
    double* x = get_frame_center(f);
    for(int d=0; d<3; d++) {
      if(x[d] > max[d]) {
	max[d] = x[d];
      }
    }
  }
  return max;
}

/*********************************
 * Tools for transform libraries *
 *********************************/

struct transform_library* create_transform_library() {
  // dynamic array for storing transforms
  struct transform_library* lib = malloc(sizeof(struct transform_library));
  assert(lib);

  lib->data = malloc(sizeof(struct transform*));
  assert(lib->data);
  lib->len = 0;
  lib->cap = 1;

  return lib;
}

void free_transform_library(struct transform_library* lib) {
  assert(lib);  
  for(int i=0; i<transform_library_size(lib); i++) {
    struct transform* T = get_transform(lib,i);
    free_transform(T);
  }
  free(lib->data);
  free(lib);
}

int transform_library_size(struct transform_library* lib) {
  assert(lib);
  return lib->len;
}

void resize_transform_library(struct transform_library* lib) {
  // double the capacity
  int new_cap = 2 * lib->cap;  
  struct transform** new_data = malloc(new_cap * sizeof(struct transform*));
  assert(new_data);
  // copy over the data to new allocation
  for(int i = 0; i < lib->len; i++) {
    new_data[i] = lib->data[i];
  }
  // free the old data, update
  free(lib->data);
  lib->data = new_data;
  lib->cap = new_cap;
}

void remove_transform(struct transform_library* lib, int idx) {
  assert(lib);
  assert((idx < lib->len && idx >= 0) || idx == -1);

  // Let users specify idx = -1 to indicate the end of the array.
  if (idx == -1) {
    idx = lib->len - 1;
  }

  // update indices after removing one at idx
  for (int i = idx; i < lib->len - 1; i++) {
    lib->data[i] = lib->data[i+1];
  }

  lib->len--;
}

struct transform* get_transform(struct transform_library* lib, int idx) {
  assert(lib);
  assert((0 <= idx && idx < lib->len) || idx == -1);

  if (idx == -1) { //  idx = -1 grabs from the end of the array.
    idx = lib->len - 1;
  }

  return lib->data[idx];
}

struct transform* get_random_transform(struct transform_library* lib) {
  int r = rand() % transform_library_size(lib);
  struct transform* T = get_transform(lib,r);
  return T;
}

void append_transform(struct transform_library* lib, struct transform* T) {
  assert(lib);
  int idx = lib->len;

  // resize if needed
  if (lib->len == lib->cap) {
    resize_transform_library(lib);
  }

  lib->data[idx] = T;
  lib->len++;
}

// combined library with both bb and bp transforms

struct library* create_library() {
  struct library* lib = malloc(sizeof(struct library));
  lib->bb = NULL;
  lib->bp = NULL;
  return lib;
}

void free_library(struct library* lib) {
  if(lib->bb) free_transform_library(lib->bb);
  if(lib->bp) free_transform_library(lib->bp);
  free(lib);
}

/**************************
 * Methods for transforms *
 **************************/

struct transform* create_transform(double val) {
  struct transform* T = malloc(sizeof(struct transform));
  T->t = create_vector(); // vector of zeros, placeholder
  T->R = create_matrix(val); // identity matrix if val=1.0
  return T;
}

void free_transform(struct transform* T) {
  assert(T);
  free_vector(T->t);
  free_matrix(T->R);
  free(T);
}

double get_transform_displacement(struct transform* T) {
  double* t = get_translation_vector(T);
  return norm(t);
}

char* get_transform_label(struct transform* T) {
  assert(T);
  return T->label;
}

void set_transform_label(struct transform* T, char* label) {
  assert(T);
  strcpy(T->label,label);  
}

double* get_translation_vector(struct transform* T) {
  return T->t;
}

void set_translation_vector(struct transform* T, double* t) {
  T->t = t;
}

double** get_rotation_matrix(struct transform* T) {
  return T->R;
}

void set_rotation_matrix(struct transform* T, double** R) {
  T->R = R;
}

struct transform* copy_transform(struct transform* T) {
  struct transform* T_copy = create_transform(1.0);
  set_vector(T_copy->t,T->t);
  set_matrix(T_copy->R,T->R);
  return T_copy;
}

struct transform* compute_transform(struct frame* f1, struct frame* f2) {
  struct transform* T = create_transform(1.0);
  double** M1_inverse = inverse(f1->M);
  double* ctr_diff = subtract(f2->ctr,f1->ctr);
  double* t = matvec_product(M1_inverse,ctr_diff);
  double** R = matmat_product(M1_inverse,f2->M);
  set_vector(T->t,t);
  set_matrix(T->R,R);
  free_vector(ctr_diff);
  free_vector(t);
  free_matrix(M1_inverse);
  free_matrix(R);
  sprintf(T->label,"%d:%d",f1->num,f2->num);
  return T;
}

struct frame* apply_transform(struct transform* T, struct frame* x1) {
  struct frame* x2 = create_frame(1.0);
  double* M_times_t = matvec_product(x1->M,T->t);
  double* ctr = add(x1->ctr,M_times_t);
  double** M = matmat_product(x1->M,T->R);
  set_vector(x2->ctr,ctr);
  set_matrix(x2->M,M);
  free_vector(M_times_t);
  free_vector(ctr);
  free_matrix(M);
  return x2;
}

struct transform* reverse_transform(struct transform* T) {
  struct transform* T_rev = create_transform(1.0);
  double** R_inv = inverse(T->R);
  double* t_neg = scale_vector(T->t,-1.0);
  double* t_rev = matvec_product(R_inv,t_neg); // transpose should work and seems faster
  set_vector(T_rev->t,t_rev);
  set_matrix(T_rev->R,R_inv);
  free_vector(t_neg);
  free_vector(t_rev);
  free_matrix(R_inv);
  return T_rev;
}
