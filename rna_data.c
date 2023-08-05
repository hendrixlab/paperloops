/* rna_data.c 
 *
 * This library contains the definitions of structures and functions for implementing
 * the RNA data structure
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "vectorfun.h"
#include "rna_data.h"

/*
 *  Data Structure Definitions 
 *  includes nodes for linked lists describing:
 *     -sequence of nucleotides (linked list of "nuc" nodes   
 *     -collection of atoms. There are three linked lists of atoms:
 *         -backbone atoms
 *         -sugar atoms
 *         -base atoms
 */

struct rna_data {
  struct rna_model** models;
  int model_count;
  int *model_locs;
  char* input_file_name;
};
  
struct rna_model {
  struct chain* first_chain;
  struct chain* last_chain;
  int chain_count;  
};

struct chain {
  struct nuc* sentinel;
  int length;
  char* seq;
  struct chain* next;
  struct chain* prev;
  char chain_id;
};

// nucleotide node
struct nuc {
  struct atom* backbone_first;
  struct atom* backbone_last;
  struct atom* sugar_first;
  struct atom* sugar_last;
  struct atom* base_first;
  struct atom* base_last;
  int num;
  char* name;
  int size;
  double bbta[6]; // backbone torsion angles
  double sta[5];  // sugar torsion angles 
  struct nuc* next;
  struct nuc* prev;
};

// atom node
struct atom {
  char* name;
  double* r;
  int num;
  char* type;
  struct atom* next;
  struct atom* prev;
  struct atom* sugar_next; // pointers just for sugar ring
  struct atom* sugar_prev;  
};

/* 
 * Node creation subroutines 
 */

struct rna_data* create_rna_data(char* input_file_name) {
  struct rna_data* RNA = malloc(sizeof(struct rna_data));
  char* alloc_string = malloc(sizeof(char)*(strlen(input_file_name)+1));
  strcpy(alloc_string,input_file_name);
  RNA->input_file_name = alloc_string;
  return RNA;
}

struct rna_model* create_model() {
  struct rna_model* RNA_model = malloc(sizeof(struct rna_model));
  RNA_model->first_chain = NULL;
  RNA_model->last_chain = NULL;
  RNA_model->chain_count = 0;
  return RNA_model;
}

struct chain* create_chain(char chain_id) {
  struct chain* RNA_chain = malloc(sizeof(struct chain));
  RNA_chain->chain_id = chain_id;
  RNA_chain->sentinel = create_nuc(0,"sentinel");
  RNA_chain->sentinel->next = RNA_chain->sentinel;
  RNA_chain->sentinel->prev = RNA_chain->sentinel;
  RNA_chain->length = 0;
  RNA_chain->next = NULL;
  RNA_chain->prev = NULL;
  return RNA_chain;
}

struct nuc* create_nuc(int nuc_num, char* nuc_name) {
  struct nuc* n = malloc(sizeof(struct nuc));
  assert(n);
  n->backbone_first = NULL;
  n->backbone_last = NULL;
  n->sugar_first = NULL;
  n->sugar_last = NULL;
  n->base_first = NULL;
  n->base_last = NULL;
  n->num = nuc_num;
  char* name_string = malloc(sizeof(char) * (strlen(nuc_name) + 1));
  strcpy(name_string,nuc_name);
  n->name = name_string;
  n->size = 0;
  n->next = NULL;
  n->prev = NULL;
  return n;
}

struct atom* create_atom(int atom_num, char* atom_name, double* r) {
  struct atom* a = malloc(sizeof(struct atom));
  assert(a);
  a->num = atom_num; 
  a->r = r;
  char* name_string = malloc(sizeof(char) * (strlen(atom_name) + 1));
  strcpy(name_string,atom_name);
  a->name = name_string;
  a->next = NULL;
  a->prev = NULL;
  a->sugar_next = NULL;
  a->sugar_prev = NULL;
  return a;
}

/* 
 * Memory management/freeing subroutines 
 */

void free_chain(struct chain* RNA_chain) {
  assert(RNA_chain);
  
  while(!chain_is_empty(RNA_chain)) {
    struct nuc* removed_nuc = get_first_nuc(RNA_chain);
    RNA_chain->sentinel = removed_nuc->next;
    removed_nuc->next->prev = RNA_chain->sentinel;
    free_nuc(removed_nuc);
  }
  free(RNA_chain);
}

void free_nuc(struct nuc* n) {
  assert(n);
  while(!nuc_backbone_is_empty(n)) {
    struct atom* removed_atom = n->backbone_first;
    n->backbone_first = removed_atom->next;
    free_atom(removed_atom);
  }
  while(!nuc_sugar_is_empty(n)) {
    struct atom* removed_atom = n->sugar_first;
    n->sugar_first = removed_atom->next;
    free_atom(removed_atom);
  }
  while(!nuc_base_is_empty(n)) {
    struct atom* removed_atom = n->base_first;
    n->base_first = removed_atom->next;
    free_atom(removed_atom);
  }  
  free(n);
}

void free_atom(struct atom* a) {
  free(a->r);
  free(a->name);
  free(a);
}

/*******************************
 * EMPTY DATA STRUCTURE CHECKS *
 *******************************/

int model_is_empty(struct rna_model* RNA_model) {
  // checks to see if the model has no chain defined, but does not check if chain is empty
  assert(RNA_model);  
  return RNA_model->first_chain == NULL;
}

int chain_is_empty(struct chain* RNA_chain) {
  // checks to see if the chain has no nucleotide defined.
  assert(RNA_chain);
  return RNA_chain->sentinel->next == RNA_chain->sentinel;
}

int nuc_backbone_is_empty(struct nuc* n) {
  assert(n);
  return n->backbone_first == NULL;
}

int nuc_sugar_is_empty(struct nuc* n) {
  assert(n);
  return n->sugar_first == NULL;
}
int nuc_base_is_empty(struct nuc* n) {
  assert(n);
  return n->base_first == NULL;
}

/* 
 * Subroutines to set/get elements to the data structure 
 */

void set_model_count(struct rna_data* RNA, int model_count) {
  RNA->model_count = model_count;
}

int get_model_count(struct rna_data* RNA) {
  return RNA->model_count;
}

void allocate_model_locs(struct rna_data* RNA) {
  RNA->model_locs = malloc(sizeof(int) * (RNA->model_count));
}

void set_model_loc(struct rna_data* RNA, int model_num, int offset) {
  RNA->model_locs[model_num] = offset;
}

int get_model_loc(struct rna_data* RNA, int model_num) {
  return RNA->model_locs[model_num];
}

void allocate_rna_models(struct rna_data* RNA) {
  RNA->models = malloc(sizeof(struct rna_model*) * (RNA->model_count));
}

void set_rna_model(struct rna_data* RNA, int model_num, struct rna_model* RNA_model) {
  RNA->models[model_num] = RNA_model;
}

int get_model_chain_count(struct rna_model* RNA_model) {
  return RNA_model->chain_count;
}

void add_chain(struct rna_model* RNA_model, struct chain* new_chain) {
  if(model_is_empty(RNA_model)) {
    RNA_model->first_chain = new_chain;  
  } else {
    new_chain->prev = RNA_model->last_chain;  
    RNA_model->last_chain->next = new_chain; // update "next" for current "last"  
  }
  RNA_model->last_chain = new_chain; // update the last node
  RNA_model->chain_count++;
}

void add_nuc(struct rna_model* RNA_model, struct nuc* n) {
  RNA_model->last_chain->sentinel->prev->next = n;
  n->prev = RNA_model->last_chain->sentinel->prev; // oddly, ->prev is the last nucleotide with circular DLL
  RNA_model->last_chain->sentinel->prev = n;
  n->next = RNA_model->last_chain->sentinel;
  RNA_model->last_chain->length++;
}

void add_atom(struct rna_model* RNA_model, struct atom* a) {
  assert(RNA_model);
  //printf("adding: %s\n",a->name);
  if(is_backbone_atom(a)) {
    // atom is part of backbone
    //printf("adding to backbone\n");
    add_backbone_atom(RNA_model,a);
  }
  if(is_sugar_atom(a)) {
    // atom is part of sugar
    //printf("adding to sugar\n");
    add_sugar_atom(RNA_model,a);
  }
  if((!is_backbone_atom(a))&&(!is_sugar_atom(a))) {
    //atom is part of the base
    //printf("adding to base\n");
    add_base_atom(RNA_model,a);
  }
}

int is_backbone_atom(struct atom* a) {
  // this function could be faster with a hash
  assert(a);
  char* backbone_atoms[] = {"P","O5'","C5'","C4'","C3'","O3'"};
  int len = 6; // hard-coded due to variable length strings
  for(int i=0; i<len; i++) {
    if(strcmp(backbone_atoms[i],a->name) == 0) {
      return 1;
    }
  }
  return 0;
}

int is_sugar_atom(struct atom* a) {
  // this function could be faster with a hash
  // not sure whether to double-add the C4' and C3' to both backbone and sugar
  char* sugar_atoms[] = {"C1'","O4'","C4'","C3'","C2'"};
  int len = 5; // hard-coded 
  for(int i=0; i<len; i++) {
    if(strcmp(sugar_atoms[i],a->name) == 0) {
      return 1;
    }
  }
  return 0;
}

void add_backbone_atom(struct rna_model* RNA_model, struct atom* a) {
  a->type = "backbone";
  if(model_is_empty(RNA_model)) {   // first if there is a nucleotide
    printf("error: attempted to add a backbone atom without a chain for the model.\n");
    exit(1);
  } else if(chain_is_empty(get_first_chain(RNA_model))) {
    printf("error: attempted to add a backbone atom without a nucleotide on the first chain.\n");
    exit(1);    
  } else {
    // next check if there is an atom attached to first nuc
    if(RNA_model->last_chain->sentinel->prev->backbone_first != NULL) {
      RNA_model->last_chain->sentinel->prev->backbone_last->next = a;
      a->prev = RNA_model->last_chain->sentinel->prev->backbone_last;
      RNA_model->last_chain->sentinel->prev->backbone_last = a;
      RNA_model->last_chain->sentinel->prev->size++;
    } else {   // this is the first atom
      RNA_model->last_chain->sentinel->prev->backbone_first = a;
      RNA_model->last_chain->sentinel->prev->backbone_last = a;
    }
  }
  RNA_model->last_chain->sentinel->prev->size++;
}

void add_sugar_atom(struct rna_model* RNA_model, struct atom* a) {
  a->type = "sugar";
  if(model_is_empty(RNA_model)) {  // check if there is a nucleotide
    printf("error: attempted to add a sugar atom without a nuc in place.\n");
    exit(1);
  } else if(chain_is_empty(get_first_chain(RNA_model))) {
    printf("error: attempted to add a backbone atom without a nucleotide on the first chain.\n");
    exit(1);    
  } else {  // next check if there is an atom attached to first nuc
    // **note: sugar atoms use separate pointers 
    if(RNA_model->last_chain->sentinel->prev->sugar_first != NULL) {
      RNA_model->last_chain->sentinel->prev->sugar_last->sugar_next = a;
      a->sugar_prev = RNA_model->last_chain->sentinel->prev->sugar_last;
      RNA_model->last_chain->sentinel->prev->sugar_last = a;
      RNA_model->last_chain->sentinel->prev->size++;
    } else {
      // this is the first atom
      RNA_model->last_chain->sentinel->prev->sugar_first = a;
      RNA_model->last_chain->sentinel->prev->sugar_last = a;
    }
  }
  RNA_model->last_chain->sentinel->prev->size++;
}

void add_base_atom(struct rna_model* RNA_model, struct atom* a) {
  a->type = "base";
  if(model_is_empty(RNA_model)) { // check if there is a nucleotide
    printf("error: attempted to add a base atom without a nuc in place.\n");
    exit(1);
  } else if(chain_is_empty(get_first_chain(RNA_model))) {
    printf("error: attempted to add a backbone atom without a nucleotide on the first chain.\n");
    exit(1);    
  } else {
    // next check if there is an atom attached to first nuc
    if(RNA_model->last_chain->sentinel->prev->base_first != NULL) {
      RNA_model->last_chain->sentinel->prev->base_last->next = a;
      a->prev = RNA_model->last_chain->sentinel->prev->base_last;
      RNA_model->last_chain->sentinel->prev->base_last = a;
      RNA_model->last_chain->sentinel->prev->size++;
    } else {
      // this is the first atom
      RNA_model->last_chain->sentinel->prev->base_first = a;
      RNA_model->last_chain->sentinel->prev->base_last = a;
    }
  }
  RNA_model->last_chain->sentinel->prev->size++;
}

/******************************************************** 
 * Functions to access elements from RNA_model data structure *
 ********************************************************/

struct chain* get_first_chain(struct rna_model* RNA_model) {
  assert(RNA_model);
  return RNA_model->first_chain;
}

//given an RNA_chain, return the next one, could be NULL, checked elsewhere
struct chain* next_chain(struct chain* RNA_chain) {
  assert(RNA_chain);
  return RNA_chain->next;
}

int new_chain_id(struct rna_model* RNA_model, char chain_id) {
  assert(RNA_model);
  if(model_is_empty(RNA_model)) {
    return 1;
  } else {
    return RNA_model->last_chain->chain_id != chain_id; // chain_id is one character long
  }
}

// this function checks if the last nucleotide in RNA_model has the number "nuc_num"
// used to avoid adding a nucleotide more than once 
int new_nuc_num(struct rna_model* RNA_model, int nuc_num) {
  assert(RNA_model);
  if(model_is_empty(RNA_model)) {
    return 1;
  } else {
    // first nuc is sentinel->prev for circular DLL
    return RNA_model->last_chain->sentinel->prev->num != nuc_num; // return true if nuc_num is a new nuc_num
  }
}

// returns a pointer to the first nucleotide (nuc) in an RNA_chain 
struct nuc* get_first_nuc(struct chain* RNA_chain) {
  assert(RNA_chain);
  if(!chain_is_empty(RNA_chain)) {
    return RNA_chain->sentinel->next;
  }
  printf("Error: No nuc found for the RNA_model\n");
  exit(1);
}

// returns a pointer to a nuc by searching for the serial number "nuc_num"
struct nuc* get_nuc_by_num(struct chain* RNA_chain, int nuc_num) {
  assert(RNA_chain);
  if(chain_is_empty(RNA_chain)) {
    struct nuc* n = RNA_chain->sentinel->next;
    while(before_end(RNA_chain,n)) {
      if(nuc_num == n->num) {
	return n;
      }
      n = next_nuc(RNA_chain,n);
    }
  }
  printf("No nuc found with number %d\n",nuc_num);
  exit(1);
}

// checks if the nuc "n" is within the RNA_model (not the sentinel)
int before_end(struct chain* RNA_chain, struct nuc* n) {
  assert(RNA_chain);
  assert(n);
  return n != RNA_chain->sentinel;
}

//given a nuc return the next nuc in the RNA_model sequence
struct nuc* next_nuc(struct chain* RNA_chain, struct nuc* n) {
  assert(n);
  return n->next;
}

//given a nuc return the next nuc in the RNA_model sequence
struct nuc* prev_nuc(struct chain* RNA_chain, struct nuc* n) {
  assert(n);
  return n->prev;
}

// returns a pointer to the previous backbone atom O3' relative to n
struct atom* prev_nuc_last_backbone_atom(struct chain* RNA_chain, struct nuc* n) {
  assert(RNA_chain);
  assert(n);
  if(!before_end(RNA_chain,prev_nuc(RNA_chain,n))) {
    return NULL;
  } else {
    return n->prev->backbone_last;
  }
}

// returns a pointer to the next backbone atom P relative to n+1
struct atom* next_nuc_first_backbone_atom(struct chain* RNA_chain, struct nuc* n) {
  assert(RNA_chain);
  assert(n);
  if(!before_end(RNA_chain,next_nuc(RNA_chain,n))) {
    return NULL;
  } else {
    return n->next->backbone_first;
  }
}

// returns a pointer to the first backbone atom for a given nucleotide
struct atom* get_first_backbone_atom(struct nuc* n) {
  assert(n);
  if(n->backbone_first) {
    struct atom* a = n->backbone_first;
    return a;    
  } else {
    printf("Error: No backbone atom found.\n");
  }
}

// returns a pointer to the last backbone atom for a given nucleotide
struct atom* get_last_backbone_atom(struct nuc* n) {
  assert(n);
  if(n->backbone_first) {
    struct atom* a = n->backbone_last;
    return a;    
  } else {
    printf("Error: No backbone atom found.\n");
  }
}


// returns a pointer to the first sugar atom
struct atom* get_first_sugar_atom(struct nuc* n) {
  assert(n);
  if(n->sugar_first) {
    struct atom* a = n->sugar_first;
    return a;    
  } else {
    printf("Error: No backbone atom found.\n");
  }
}

// returns a pointer to the next atom in a set given pointer to atom
// e.g. the next backbone atom
struct atom* next_atom(struct atom* a) {
  assert(a);
  return a->next;
}

// returns a pointer to the next atom in a set given pointer to atom
// e.g. the next backbone atom
struct atom* prev_atom(struct atom* a) {
  assert(a);
  return a->prev;
}

// returns a pointer to the next atom in a set given pointer to atom
// e.g. the next backbone atom
struct atom* next_sugar_atom(struct atom* a) {
  assert(a);
  return a->sugar_next;
}

// returns a pointer to the next atom in a set given pointer to atom
// e.g. the next backbone atom
struct atom* prev_sugar_atom(struct atom* a) {
  assert(a);
  return a->sugar_prev;
}

// returns a pointer to an atom by searching for the serial number "atom_num"
struct atom* get_atom_by_num(struct nuc* n, int atom_num) {
  assert(n);
  if(n->backbone_first) {
    struct atom* a = n->backbone_first;
    while(a) {
      if(atom_num == a->num) {
	return a;
      }
      a = a->next;
    }
  }
  if(n->sugar_first) {
    struct atom* a = n->sugar_first;
    while(a) {
      if(atom_num == a->num) {
	return a;
      }
      a = a->sugar_next;
    }
  }
  if(n->base_first) {
    struct atom* a = n->base_first;
    while(a) {
      if(atom_num == a->num) {
	return a;
      }
      a = a->next;
    }
  }
  printf("No atom found with number %d\n",atom_num);
  exit(1);
}

// prints all sugar atoms for a nucleotide
void print_sugar_atoms(struct nuc* n) {
  assert(n);
  printf("(%s)->",n->name);
  if(n->sugar_first) {
    struct atom* a = n->sugar_first;
    int count = 0;
    while(count<5) {
      printf("[%s]",a->name);
      a = a->sugar_next;
      count++;
    }
  } else {
    printf("-|");
  }
  printf("\n");
}


// returns a pointer to an atom by searching by name
struct atom* get_sugar_atom_by_name(struct nuc* n, char* atom_name) {
  assert(n);
  if(n->sugar_first) {
    struct atom* a = n->sugar_first;
    while(a) {
      if(strcmp(atom_name,a->name) == 0) {
	return a;
      }
      a = a->sugar_next;
    }
  }
  printf("No atom found with name: %s\n",atom_name);
  exit(1);
}

// returns a pointer to an atom by searching by name
struct atom* remove_sugar_atom_by_name(struct nuc* n, char* atom_name) {
  assert(n);
  if(n->sugar_first) {
    struct atom* a = n->sugar_first;
    while(a) {
      if(strcmp(atom_name,a->name) == 0) {
	if(a == n->sugar_first) {
	  if(a == n->sugar_last) { // both first and last -> only one
	    n->sugar_first = NULL;
	    n->sugar_last = NULL;
	  } else { // at least one other node after
	    n->sugar_first = a->sugar_next;
	    n->sugar_first->prev = NULL;
	  }
	} else if(a == n->sugar_last) {  // last but not first
	  a->sugar_prev->sugar_next = NULL;
	  n->sugar_last = a->sugar_prev;
	} else {
	  a->sugar_prev->sugar_next = a->sugar_next;
	  a->sugar_next->sugar_prev = a->sugar_prev;
	}
	a->sugar_next = NULL;
	a->sugar_prev = NULL;
	return a;
      }
      a = a->sugar_next;
    }
  }
  return NULL;
}

void sort_sugar_atoms(struct nuc* n) {
  struct atom* new_sugar_first = NULL;
  struct atom* new_sugar_last = NULL;
  //print_sugar_atoms(n);
  char* sugar_atoms[] = {"C1'","O4'","C4'","C3'","C2'"};
  int len = 5; // hard-coded
  for(int i=0; i<len; i++) {
    struct atom* a = remove_sugar_atom_by_name(n,sugar_atoms[i]);
    if(!a) {
      // if atom is missing, create empty placeholder
      a = create_atom(0,sugar_atoms[i],NULL);
    }
    if(i==0) {
      a->sugar_next = a;
      a->sugar_prev = a;
      new_sugar_first = a;
      new_sugar_last = a;
    } else {
      a->sugar_prev = new_sugar_last;
      a->sugar_next = new_sugar_first;
      new_sugar_last->sugar_next = a;
      new_sugar_last = a;
      new_sugar_first->sugar_prev = new_sugar_last;
    }
  }
  n->sugar_first = new_sugar_first;
  n->sugar_last = new_sugar_last;
  //printf("New sort:");
  //print_sugar_atoms(n);
}

// retrieve the serial number of a nucletide 
int nuc_num(struct nuc* n) {
  assert(n);
  return n->num;
}

// retrieve the name of a nucleotide 
char* nuc_name(struct nuc* n) {
  assert(n);
  return n->name;
}

// retrieve the serial number of an atom
int atom_num(struct atom* a) {
  assert(a);
  return a->num;
}

// retrieve the name of an atom
char* atom_name(struct atom* a) {
  assert(a);
  return a->name;
}
 
// retrieve the type of an atom
char* atom_type(struct atom* a) {
  assert(a);
  return a->type;
}

// retrieve the atom coordinate
double* atom_coord(struct atom* a) {
  assert(a);
  return a->r;
}

/******************************
 * Torsion angle calculations *
 ******************************/

double torsion1(struct atom* p1, struct atom* p2, struct atom* p3, struct atom *p4) {
  // method 1: most intutive method, from: https://x3dna.org/highlights/how-to-calculate-torsion-angle
  if(p1 && p2 && p3 && p4) {
    if(p1->r && p2->r && p3->r && p4->r) {
      double* vec_a = malloc(3 * sizeof(double));
      double* vec_b = malloc(3 * sizeof(double));
      double* vec_c = malloc(3 * sizeof(double));
      vec_a = normalize(subtract(p1->r,p2->r)); // note p1-p2 so it is from p2 to p1
      vec_b = normalize(subtract(p3->r,p2->r));
      vec_c = normalize(subtract(p4->r,p3->r));
      double* proj_a = orthonormal_projection(vec_a,vec_b);
      double* proj_c = orthonormal_projection(vec_c,vec_b);
      double angle = get_angle(proj_a,proj_c);
      double sign = dot(cross(proj_a,proj_c),vec_b);
      if(sign < 0) {
	angle = -angle;
      }
      return angle;
    } else {
      return NAN;
    }
  } else {
    return NAN;
  }
}

double torsion2(struct atom* p1, struct atom* p2, struct atom* p3, struct atom *p4) {
  // method 2: simple with cross products, described in: https://www.math.fsu.edu/~quine/MB_10/6_torsion.pdf
  // I added a sign correction based on cross product with vec_b --> same result as method 1
  if(p1 && p2 && p3 && p4) {
    double* vec_a = malloc(3 * sizeof(double));
    double* vec_b = malloc(3 * sizeof(double));
    double* vec_c = malloc(3 * sizeof(double));
    vec_a = normalize(subtract(p2->r,p1->r)); // note p1-p2 so it is from p2 to p1
    vec_b = normalize(subtract(p3->r,p2->r));
    vec_c = normalize(subtract(p4->r,p3->r));
    double angle = get_angle(cross(vec_a,vec_b),cross(vec_b,vec_c));
    double sign = dot(cross(cross(vec_a,vec_b),cross(vec_b,vec_c)),vec_b);
    if(sign < 0) {
      angle = -angle;
    }
    return angle;
  } else {
    return NAN;
  }
}

double torsion3(struct atom* p1, struct atom* p2, struct atom* p3, struct atom *p4) {
  // method 3: https://www.math.fsu.edu/~quine/MB_10/6_torsion.pdf
  // seemed off by an additive factor of pi
  if(p1 && p2 && p3 && p4) {
    double* vec_a = malloc(3 * sizeof(double));
    double* vec_b = malloc(3 * sizeof(double));
    double* vec_c = malloc(3 * sizeof(double));
    vec_a = normalize(subtract(p2->r,p1->r)); // note p1-p2 so it is from p2 to p1
    vec_b = normalize(subtract(p3->r,p2->r));
    vec_c = normalize(subtract(p4->r,p3->r));  
    //double angle = argument(dot(cross(vec_b,vec_c),cross(vec_a,vec_b)),dot(cross(vec_b,vec_c),cross(vec_b,cross(vec_b,vec_c)))/norm(vec_b));
    double angle = argument(dot(vec_a,vec_b)*dot(vec_b,vec_c)-dot(vec_b,vec_b)*dot(vec_a,vec_c),norm(vec_b)*dot(vec_a,cross(vec_b,vec_c)));
    return angle;
  } else {
    return NAN;
  }
}

void compute_backbone_torsions(struct rna_data* RNA) {
  for(int i=0; i<RNA->model_count; i++) {
    compute_model_backbone_torsions(RNA->models[i]);
  }
}

void compute_sugar_torsions(struct rna_data* RNA) {
  for(int i=0; i<RNA->model_count; i++) {
    compute_model_sugar_torsions(RNA->models[i]);
  }
}

void compute_model_backbone_torsions(struct rna_model* RNA_model) {
  struct chain* RNA_chain = get_first_chain(RNA_model);
  while(RNA_chain) {
    struct nuc* n = get_first_nuc(RNA_chain);
    int FIRST = 1;
    /* potential bug: should this be after_end only? --> may explain 0.0 angles*/
    while(before_end(RNA_chain,n)) { 
      struct atom* p1 = prev_nuc_last_backbone_atom(RNA_chain,n);
      struct atom* p2 = get_first_backbone_atom(n);
      struct atom* p3 = next_atom(p2);
      struct atom* p4 = next_atom(p3);
      int count;
      if(strcmp(atom_name(p2),"P") == 0) {
	count = 0;
      } else if(strcmp(atom_name(p2),"O5'") == 0) {
	n->bbta[0] = NAN;
	count = 1;
      } else {
	printf("Unexpected first atom in backbone: \"%s\"\n",atom_name(p2));
	exit(1);
      }
      int CONTINUE = 1;
      int NEXT_NUC = 0;
      while(CONTINUE) {
	n->bbta[count] = torsion1(p1,p2,p3,p4); // store torsion angle
	p1 = p2; // --> next quadruple
	p2 = p3;
	p3 = p4;
	if(p4 && next_atom(p4)) {
	  p4 = next_atom(p4);
	  if(NEXT_NUC > 1) {
	    CONTINUE = 0; // allow two atoms from next nucleotide P & O5'
	  } else if(NEXT_NUC > 0) {
	    NEXT_NUC++;
	  }
	} else {
	  if(before_end(RNA_chain,next_nuc(RNA_chain,n))) {
	    p4 = next_nuc_first_backbone_atom(RNA_chain,n);
	    NEXT_NUC++;
	  } else {
	    CONTINUE = 0;
	  }
	}
	count++;
      }
      n = next_nuc(RNA_chain,n);
    }
    RNA_chain = next_chain(RNA_chain);
  }
}

void compute_model_sugar_torsions(struct rna_model* RNA_model) {
  struct chain* RNA_chain = get_first_chain(RNA_model);
  while(RNA_chain) {
    struct nuc* n = get_first_nuc(RNA_chain);  
    while(before_end(RNA_chain,n)) {
      // sort atoms in bonded order
      sort_sugar_atoms(n);
      int count = 0;
      struct atom* p2 = get_first_sugar_atom(n);
      struct atom* p1 = prev_sugar_atom(p2);
      struct atom* p3 = next_sugar_atom(p2);
      struct atom* p4 = next_sugar_atom(p3);;
      while(count<5) {
	n->sta[count] = torsion1(p1,p2,p3,p4); // store torsion angle    
	p1 = p2;
	p2 = p3;
	p3 = p4;
	p4 = next_sugar_atom(p3);
	count++;
      }
      n = next_nuc(RNA_chain,n);
    }
    RNA_chain = next_chain(RNA_chain);
  }
}

void print_backbone_torsions(struct rna_data* RNA) {
  char ext[] = ".backbone.txt";
  char output_file[100];
  strcpy(output_file,RNA->input_file_name);
  strcat(output_file,ext);
  FILE* OUTPUT = fopen(output_file,"w");  

  char* angle_name[] = {"alpha","beta","gamma","delta","epsilon","zeta"};
  /* header */
  fprintf(OUTPUT,"model\tchain\tpos\tnt");
  for(int j=0; j<6; j++) {
    fprintf(OUTPUT,"\t%s",angle_name[j]);
  }
  fprintf(OUTPUT,"\n");
  /* data table */  
  for(int i=0; i<RNA->model_count; i++) {
    struct chain* RNA_chain = get_first_chain(RNA->models[i]);
    while(RNA_chain) {
      struct nuc* n = get_first_nuc(RNA_chain);
      int FIRST = 1;
      while(before_end(RNA_chain,n)) {
	fprintf(OUTPUT,"%d\t%c\t%d\t%s",i,RNA_chain->chain_id,n->num,n->name);
	for(int j=0; j<6; j++) {
	  fprintf(OUTPUT,"\t%.2f",n->bbta[j]);
	}
	fprintf(OUTPUT,"\n");
	n = next_nuc(RNA_chain,n);
      }
      RNA_chain = next_chain(RNA_chain);
    }
  }
  fclose(OUTPUT);
}

void print_sugar_torsions(struct rna_data* RNA) {  
  char ext[] = ".sugar.txt";
  char output_file[100];
  strcpy(output_file,RNA->input_file_name);
  strcat(output_file,ext);
  FILE* OUTPUT = fopen(output_file,"w");  

  char* angle_name[] = {"v1","v2","v3","v4","v5"};
  int len = 5;
  /* header */
  fprintf(OUTPUT,"model\tchain\tpos\tnt");
  for(int j=0; j<len; j++) {
    fprintf(OUTPUT,"\t%s",angle_name[j]);
  }
  fprintf(OUTPUT,"\n");
  /* data table */
  for(int i=0; i<RNA->model_count; i++) {
    struct chain* RNA_chain = get_first_chain(RNA->models[i]);
    while(RNA_chain) {
      struct nuc* n = get_first_nuc(RNA_chain);      
      int FIRST = 1;  
      while(before_end(RNA_chain,n)) {
	fprintf(OUTPUT,"%d\t%d\t%s",i,n->num,n->name);
	for(int j=0; j<len; j++) {
	  fprintf(OUTPUT,"\t%.2f",n->sta[j]);
	}
	fprintf(OUTPUT,"\n");
	n = next_nuc(RNA_chain,n);
      }
      RNA_chain = next_chain(RNA_chain);
    }
  }
  fclose(OUTPUT);  
}

/* Functions for extracting the RNA_model sequence */

char* get_rna_sequence(struct rna_data* RNA) {
  assert(RNA);
  if(RNA->model_count > 0) {
    for(int i=0; i<RNA->model_count; i++) {      
      if(RNA->models[i]->chain_count > 0) {
	struct chain* RNA_chain = get_first_chain(RNA->models[i]);
	while(RNA_chain) {
	  char* seq = get_chain_sequence(RNA_chain);
	  printf(">%c\n%s\n",RNA_chain->chain_id,seq);
	  RNA_chain = next_chain(RNA_chain);
	}
      }
    }
  }
}

char* set_chain_sequence(struct chain* RNA_chain) {
  assert(RNA_chain);
  char* sequence; 
  if(RNA_chain->length != 0) {
    sequence = malloc(sizeof(char) * RNA_chain->length);
    struct nuc* n = RNA_chain->sentinel->next;
    strcpy(sequence,"");
    while(before_end(RNA_chain,n)) {
      strcat(sequence,n->name);
      n = next_nuc(RNA_chain,n);
    }
  } else {
    printf("Error: RNA_chain has length zero.\n");
    exit(0);
  }
  RNA_chain->seq = sequence;
}

char* get_chain_sequence(struct chain* RNA_chain) {
  assert(RNA_chain);
  if(RNA_chain->seq == NULL) {
    set_chain_sequence(RNA_chain);
  } else if(strlen(RNA_chain->seq) != RNA_chain->length) {
    set_chain_sequence(RNA_chain);    
  }
  return RNA_chain->seq;
}
