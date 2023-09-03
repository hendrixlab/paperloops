 /*
 * This library contains the definitions of structures and functions for transforms between 
 * adjacent nucleotides
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "transform.h"
#include "rna_data.h"
#include "vectorfun.h"

/*
 * data structures for the transform and for the library of transforms
 */

struct frame {
  double* ctr;
  double** M;
};

struct transform {
  double* t;
  double** R;
};

struct library {
  struct transform** data;
  int len;
  int cap;
};

/*
 * Tools for building/adding/removing elements from library
 */

struct library* create_library() {
  struct library* transform_set = malloc(sizeof(struct library));
  assert(transform_set);

  transform_set->data = malloc(sizeof(struct transform*));
  assert(transform_set->data);
  transform_set->len = 0;
  transform_set->cap = 1;

  return transform_set;
}

void free_library(struct library* transform_set) {
  assert(transform_set);
  free(transform_set->data);
  free(transform_set);
}

int library_size(struct library* transform_set) {
  assert(transform_set);
  return transform_set->len;
}

void resize_library(struct library* transform_set) {
  // double the capacity
  int new_cap = 2 * transform_set->cap;  
  struct transform** new_data = malloc(new_cap * sizeof(struct transform*));
  assert(new_data);
  // copy over the data to new allocation
  for (int i = 0; i < transform_set->len; i++) {
    new_data[i] = transform_set->data[i];
  }
  // free the old data, update
  free(transform_set->data);
  transform_set->data = new_data;
  transform_set->cap = new_cap;
}

void remove_transform(struct library* transform_set, int idx) {
  assert(transform_set);
  assert((idx < transform_set->len && idx >= 0) || idx == -1);

  // Let users specify idx = -1 to indicate the end of the array.
  if (idx == -1) {
    idx = transform_set->len - 1;
  }

  // update indices after removing one at idx
  for (int i = idx; i < transform_set->len - 1; i++) {
    transform_set->data[i] = transform_set->data[i+1];
  }

  transform_set->len--;
}

struct transform* get_transform(struct library* transform_set, int idx) {
  assert(transform_set);
  assert((idx < transform_set->len && idx >= 0) || idx == -1);

  // Let users specify idx = -1 to indicate the end of the array.
  if (idx == -1) {
    idx = transform_set->len - 1;
  }

  return transform_set->data[idx];
}

void add_transform(struct library* transform_set, struct transform* T) {
  assert(transform_set);
  int idx = transform_set->len;

  // resize if needed
  if (transform_set->len == transform_set->cap) {
    resize_library(transform_set);
  }

  transform_set->data[idx] = T;
  transform_set->len++;
}

/*
 * Methods for building the transform library
 */

struct frame* create_frame() {
  struct frame* f = malloc(sizeof(struct frame*));
  f->ctr = malloc(3*sizeof(double));
  f->M = malloc(3*sizeof(double*));
  f->M[0] = malloc(3*sizeof(double));
  f->M[1] = malloc(3*sizeof(double));
  f->M[2] = malloc(3*sizeof(double));
  return f;
}

struct frame* get_frame(struct nuc* n) {
  struct frame* f = create_frame();
  struct atom* a = get_backbone_atom_by_name(n,"C5'"); 
  struct atom* b = get_backbone_atom_by_name(n,"C4'"); 
  struct atom* c = get_backbone_atom_by_name(n,"C3'");
  double* z = normalize(subtract(atom_coord(b),atom_coord(a)));
  double* y = normalize(cross(z,subtract(atom_coord(c),atom_coord(b))));
  double* x = normalize(cross(y,z));
  f->ctr = atom_coord(a);
  f->M[0] = x;
  f->M[1] = y;
  f->M[2] = z;
  return f;
}

struct transform* compute_transform(struct frame* f1, struct frame* f2) {
  struct transform* T = malloc(sizeof(struct transform));
  T->t = malloc(3*sizeof(double));
  T->R = malloc(3*sizeof(double*));
  T->R[0] = malloc(3*sizeof(double));
  T->R[1] = malloc(3*sizeof(double));
  T->R[2] = malloc(3*sizeof(double));
  double** M1_inverse = inverse(f1->M);
  T->t = matvec_product(M1_inverse,subtract(f2->ctr,f1->ctr));
  T->R = matmat_product(M1_inverse,f2->M);
  return T;
}

int valid_dinucleotide(struct nuc* n1, struct nuc* n2) {
  struct atom* a1 = get_last_backbone_atom(n1);
  struct atom* a2 = get_first_backbone_atom(n2);
  if(a1 && a2) {
    if((strcmp(atom_name(a1),"O3'") == 0)&&(strcmp(atom_name(a2),"P") == 0)) {
      double d = norm(subtract(atom_coord(a1),atom_coord(a2)));
      if(d < 2) {
	return 1;
      }    
    } else {
      printf("Unexpected atoms found: %s %s\n",atom_name(a1),atom_name(a2));
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
      while(before_end(chain,n2)) {
	if(valid_dinucleotide(n1,n2)) {
	  //printf("found pair %d: %d %d (%s,%s)\n",count,nuc_num(n1),nuc_num(n2),nuc_name(n1),nuc_name(n2));
	  struct frame* f1 = get_frame(n1);
	  struct frame* f2 = get_frame(n2);
	  struct transform* T = compute_transform(f1,f2);
	  add_transform(backbone_library,T);
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
      //set_chain_sequence(chain1);
      //printf("chain1: %s\n",get_chain_sequence(chain1));
      struct nuc* n1 = get_nuc_by_idx(chain1,pos1);
      struct nuc* n2 = get_nuc_by_idx(chain2,pos2);
      struct frame* f1 = get_frame(n1);
      struct frame* f2 = get_frame(n2);
      struct transform* T = compute_transform(f1,f2);
      add_transform(basepair_library,T);
    }
  }
  printf("found %d base pairs.\n",library_size(basepair_library));
  return basepair_library;
}


