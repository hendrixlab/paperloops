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
#include "motifs.h"
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

/*
 * Methods for building the transform library
 */

struct frame* create_frame() {
  struct frame* f = malloc(sizeof(struct frame*));
  f->ctr = create_vector(); // vector of zeros, placeholer
  f->M = create_matrix(); // identity matrix, placeholder
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

struct transform* create_transform() {
  struct transform* T = malloc(sizeof(struct transform));
  T->t = create_vector(); // vector of zeros, placeholder
  T->R = create_matrix(); // identity matrix, placeholder
  return T;
}

struct transform* compute_transform(struct frame* f1, struct frame* f2) {
  struct transform* T = create_transform();
  double** M1_inverse = inverse(f1->M);
  T->t = matvec_product(M1_inverse,subtract(f2->ctr,f1->ctr));
  T->R = matmat_product(M1_inverse,f2->M);
  return T;
}

struct transform* apply_transform(struct transform* x1, struct transform* T) {
  struct transform* x2 = create_transform();
  x2->t = add(x1->t,matvec_product(x1->R,T->t));
  x2->R = matmat_product(x1->R,T->R);
  return x2;
}

struct transform* reverse_transform(struct transform* T) {
  struct transform* T_rev = create_transform();
  T_rev->t = matvec_product(inverse(T->R),negate(T->t));;
  T_rev->R = inverse(T->R);
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
      while(before_end(chain,n2)) { // while n2 is not a sentinel & is a valid nt
	if(valid_dinucleotide(n1,n2)) { // checks that no gap b/w n1 and n2
	  struct frame* f1 = get_frame(n1);
	  struct frame* f2 = get_frame(n2);
	  struct transform* T = compute_transform(f1,f2);
	  append_transform(backbone_library,T);
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
    }
  }
  printf("found %d base pairs.\n",library_size(basepair_library));
  return basepair_library;
}

/***********************
 * sampling/simulation *
 ***********************/

struct library* get_random_trace(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib) {
  srand(1337);
  struct library* trace = create_library(); // library to store the trace
  struct transform* x = create_transform(); // the initial point

  for(int i=0; i<steptypes_length(steps); i++) {
    char* step = steptypes_get(steps,i);
    struct transform* T;
    if(strcmp(step,"bb") == 0) {
      T = get_random_transform(bb_lib);
    } else if(strcmp(step,"bp") == 0) {
      T = get_random_transform(bp_lib);
    } else {
      printf("Unknown step found!\n");
      exit(1);
    }
    x = apply_transform(x,T);
    append_transform(trace,x);
  }
  return trace;
}

struct library* get_pts_forward(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size) {
  struct library* pts_fwd = create_library(); // create the library to store the points
  for(int i=0; i<size; i++) {
    struct library* trace = get_random_trace(steps,bb_lib,bp_lib);
    struct transform* x = get_transform(trace,-1); // store the last point
    append_transform(pts_fwd,x);
  }
  return pts_fwd;
}

struct library* get_pts_reverse(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size) {
  struct library* pts_fwd = get_pts_forward(steps,bb_lib,bp_lib,size);
  struct library* pts_rev = create_library();
  for(int i=0; i<library_size(pts_fwd); i++) {
    struct transform* T = get_transform(pts_fwd,i);
    struct transform* T_rev = reverse_transform(T);
    append_transform(pts_rev,T_rev);
  }
}



