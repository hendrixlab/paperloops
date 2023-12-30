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

/**********************************************************
 * Data structures for transforms and transform libraries *
 **********************************************************/

struct frame {
  double* ctr;
  double** M;
};

struct transform {
  double* t; // translation
  double** R; // rotation
  double* EV; // Euler vector
};

struct library {
  struct transform** data;
  int len;
  int cap;
};

/*********************************
 * Tools for transform libraries *
 *********************************/

struct library* create_library() {
  struct library* lib = malloc(sizeof(struct library));
  assert(lib);

  lib->data = malloc(sizeof(struct transform*));
  assert(lib->data);
  lib->len = 0;
  lib->cap = 1;

  return lib;
}

void free_library(struct library* lib) {
  assert(lib);  
  for(int i=0; i<library_size(lib); i++) {
    struct transform* T = get_transform(lib,i);
    free_transform(T);
  }
  free(lib->data);
  free(lib);
}

int library_size(struct library* lib) {
  assert(lib);
  return lib->len;
}

void resize_library(struct library* lib) {
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

void remove_transform(struct library* lib, int idx) {
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

struct transform* get_transform(struct library* lib, int idx) {
  assert(lib);
  assert((0 <= idx && idx < lib->len) || idx == -1);

  if (idx == -1) { //  idx = -1 grabs from the end of the array.
    idx = lib->len - 1;
  }

  return lib->data[idx];
}

struct transform* get_random_transform(struct library* lib) {
  int r = rand() % library_size(lib);
  struct transform* T = get_transform(lib,r);
  return T;
}

void append_transform(struct library* lib, struct transform* T) {
  assert(lib);
  int idx = lib->len;

  // resize if needed
  if (lib->len == lib->cap) {
    resize_library(lib);
  }

  lib->data[idx] = T;
  lib->len++;
}

/**********************
 * Methods for frames *
 **********************/

struct frame* create_frame() {
  struct frame* f = malloc(sizeof(struct frame));
  f->ctr = create_vector(); // vector of zeros, placeholer
  f->M = create_matrix(1.0); // identity matrix, placeholder
  return f;
}

void free_frame(struct frame* f) {
  assert(f);
  free_vector(f->ctr);
  free_matrix(f->M);
  free(f);
}
  
struct frame* get_frame(struct nuc* n) {
  struct frame* f = create_frame();
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

/**************************
 * Methods for transforms *
 **************************/

struct transform* create_transform() {
  struct transform* T = malloc(sizeof(struct transform));
  T->t = create_vector(); // vector of zeros, placeholder
  T->R = create_matrix(1.0); // identity matrix, placeholder
  T->EV = create_vector(); // Euler vector, initially matrix of zeros
  return T;
}

void free_transform(struct transform* T) {
  assert(T);
  free_vector(T->t);
  free_matrix(T->R);
  free_vector(T->EV);
  free(T);
}

double** get_rotation_matrix(struct transform* T) {
  return T->R;
}

struct transform* copy_transform(struct transform* T) {
  struct transform* T_copy = create_transform();
  set_vector(T_copy->t,T->t);
  set_matrix(T_copy->R,T->R);
  set_vector(T_copy->EV,T->EV);
  return T_copy;
}

struct transform* compute_transform(struct frame* f1, struct frame* f2) {
  struct transform* T = create_transform();
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
  return T;
}

struct transform* apply_transform(struct transform* x1, struct transform* T) {
  struct transform* x2 = create_transform();
  double* R_times_t = matvec_product(x1->R,T->t);
  double* t = add(x1->t,R_times_t);
  double** R = matmat_product(x1->R,T->R);
  set_vector(x2->t,t);
  set_matrix(x2->R,R);
  free_vector(R_times_t);
  free_vector(t);
  free_matrix(R);
  return x2;
}

struct transform* reverse_transform(struct transform* T) {
  struct transform* T_rev = create_transform();
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

void set_euler_vector(struct transform* T) {
  // calling function in vectorfun.c
  double* EV = compute_euler_vector(T->R);
  set_vector(T->EV,EV);
  free_vector(EV);
}

double* get_euler_vector(struct transform* T) {
  return T->EV;
}

double get_T6_value(struct transform* T, int i) {
  if(0 <= i && i < 3) {
    return T->t[i];
  } else if(3 <= i && i < 6) { // hard-coded b/c should never be other than 6D
    return T->EV[i-3];
  } else {
    printf("Unexpected value in get_T6_value: %d",i);
    exit(0);
  }
}

double* get_T6_std_dev(struct library* pts) {
  int d = 6; // assuming 6D coordinates
  double* sum = malloc(sizeof(double)*d);
  double* sum_squares = malloc(sizeof(double)*d);
  double* std_dev = malloc(sizeof(double)*d);  
  int n = library_size(pts);
  
  for(int j=0; j<d; j++) {
    sum[j] = 0.0;
    sum_squares[j] = 0.0;
  }
  for(int i=0; i<n; i++) {
    struct transform* T = get_transform(pts,i);
    for(int j=0; j<d; j++) {
      double val = get_T6_value(T,j);
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

void print_transform_T6(struct transform* x) {
  int d=6;
  printf("[");
  for(int j=0; j<d; j++) {
    double x_val = get_T6_value(x,j);
    printf(" %f ",x_val);
  }
  printf("]\n");
}

/*****************************************************
 * processing biophysical data to compute transforms *
 *****************************************************/

int valid_dinucleotide(struct nuc* n1, struct nuc* n2) {
  struct atom* a1 = get_last_backbone_atom(n1);
  struct atom* a2 = get_first_backbone_atom(n2);
  if(a1 && a2) {
    if((strcmp(atom_name(a1),"O3'") == 0)&&(strcmp(atom_name(a2),"P") == 0)) {
      double* diff = subtract(atom_coord(a1),atom_coord(a2));
      double d = norm(diff);
      free_vector(diff);
      if(d < 2) {
	return 1;
      }    
    } else {
      printf("Unexpected atoms found for nucs %d and %d: %s %s\n",nuc_num(n1),nuc_num(n2),atom_name(a1),atom_name(a2));
      exit(1);
    }
  }
  return 0;
}

struct library* build_backbone_transform_library(struct rna_data* RNA) {
  assert(RNA);
  if(get_model_count(RNA) == 0) {
    printf("No models found for RNA.\n");
    exit(1);
  }
  struct library* backbone_library = create_library();
  for(int i=0; i<get_model_count(RNA); i++) {
    struct rna_model* RNA_model = get_rna_model(RNA,i); 
    struct chain* chain = get_first_chain(RNA_model);
    while(chain) {
      struct nuc* n1 = get_first_nuc(chain);
      struct nuc* n2 = next_nuc(chain,n1);
      while(before_end(chain,n2)) { // while n2 is not a sentinel & is a valid nt
	if(valid_dinucleotide(n1,n2)) { // checks that no gap b/w n1 and n2
	  struct frame* f1 = get_frame(n1);
	  struct frame* f2 = get_frame(n2);
	  struct transform* T = compute_transform(f1,f2);
	  append_transform(backbone_library,T);
	  free_frame(f1); // free memory of each frame after use
	  free_frame(f2);
	}
	n1 = n2;
	n2 = next_nuc(chain,n1);	   
      }
      chain = next_chain(chain);
    }
  }
  printf("found %d dinucleotide pairs.\n",library_size(backbone_library));
  return backbone_library;
}

struct library* build_basepair_transform_library(struct rna_data* RNA, struct bp_array* bp_array) {
  assert(RNA);
  if(get_model_count(RNA) == 0) {
    printf("No models found for RNA.\n");
    exit(1);
  }
  struct library* basepair_library = create_library();
  for(int i=0; i<get_model_count(RNA); i++) {      
    struct rna_model* RNA_model = get_rna_model(RNA,i);
    for(int j=0; j<get_bp_array_size(bp_array);j++) {
      struct bp* bp = get_bp(bp_array,j);
      struct chain* chain1 = get_chain_by_id(RNA_model,get_bp_chain_id(bp,1));
      struct chain* chain2 = get_chain_by_id(RNA_model,get_bp_chain_id(bp,2));
      int pos1 = get_bp_pos(bp,1);
      int pos2 = get_bp_pos(bp,2);
      struct nuc* n1 = get_nuc_by_idx(chain1,pos1);
      struct nuc* n2 = get_nuc_by_idx(chain2,pos2);
      struct frame* f1 = get_frame(n1);
      struct frame* f2 = get_frame(n2);
      struct transform* T = compute_transform(f1,f2);
      append_transform(basepair_library,T);
      free_frame(f1);
      free_frame(f2);
    }
  }
  printf("found %d base pairs.\n",library_size(basepair_library));
  return basepair_library;
}
