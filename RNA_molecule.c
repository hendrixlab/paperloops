/*
 * This file contains the definitions of structures and functions implementing
 * for the RNA data structure
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "RNA_molecule.h"

/*
 * This is the definition of the RNA linked list structure.  
 */

struct rna {
  struct nuc* first;
  struct nuc* last;
  int length;
};

/*
 * This is the definition of the "nuc" node for the linked list.
 */

struct nuc {
  struct atom* first;
  struct atom* last;
  int num;
  char* name;
  int size;
  struct nuc* next;
};

/*
 * This is the definition of the "atom" node for each nucleotide "nuc"
 */

struct atom {
  char* name;
  double* r;
  int num;
  struct atom* next;
};

/* creation subroutines */

struct rna* create_rna() {
  struct rna* RNA = malloc(sizeof(struct rna));
  assert(RNA);
  RNA->first = NULL;
  RNA->last = NULL;
  RNA->length = 0;
  return RNA;
}

struct nuc* create_nuc(int nuc_num, char* nuc_name) {
  struct nuc* n = malloc(sizeof(struct nuc));
  assert(n);
  n->first = NULL;
  n->last = NULL;
  n->num = nuc_num;
  char* name_string = malloc(sizeof(char) * (strlen(nuc_name) + 1));
  strcpy(name_string,nuc_name);
  n->name = name_string;
  n->size = 0;
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
  return a;
}

/* free memory subroutines */

void free_rna(struct rna* RNA) {
  assert(RNA);
  
  while(!rna_is_empty(RNA)) {
    struct nuc* removed_nuc = RNA->first;
    RNA->first = removed_nuc->next;
    free_nuc(removed_nuc);
  }
  free(RNA);
}

void free_nuc(struct nuc* n) {
  assert(n);

  while(!nuc_is_empty(n)) {
    struct atom* removed_atom = n->first;
    n->first = removed_atom->next;
    free_atom(removed_atom);
  }
  free(n);
}

void free_atom(struct atom* a) {
  free(a->r);
  free(a->name);
  free(a);
}

int rna_is_empty(struct rna* RNA) {
  assert(RNA);
  return RNA->first == NULL;
}

int nuc_is_empty(struct nuc* n) {
  assert(n);
  return n->first == NULL;
}

/* adding elements to the data structure */

void add_nuc(struct rna* RNA, struct nuc* n) {
  if(RNA->first != NULL) {
    RNA->last->next = n;
    RNA->last = n;
  } else {
    RNA->first = n;
    RNA->last = n;
  }
  RNA->length++;    
}

void add_atom(struct rna* RNA, struct atom* a) {
  // first check if there is a nucleotide
  if(RNA->first != NULL) {
    // next check if there is an atom attached to first nuc
    if(RNA->last->first != NULL) {
      RNA->last->last->next = a;
      RNA->last->last = a;
      RNA->last->size++;
    } else {
      // this is the first atom
      RNA->last->first = a;
      RNA->last->last = a;
    }
  } else {
    // error not expecting to add atom without nuc already in place
    printf("error: attempted to add an atom without a nuc in place.\n");
    exit(1);
  }
  RNA->last->size++;
}

/* functions to access elements from RNA data structure */

// this function checks if the last nucleotide in RNA has the number "nuc_num"
// used to avoid adding a nucleotide more than once 
int last_nuc_has_num(struct rna* RNA, int nuc_num) {
  assert(RNA);
  if(RNA->last != NULL) {
    return RNA->last->num == nuc_num;
  } else {
    return 0;
  }
}

// returns a pointer to a nuc by searching for the serial number "nuc_num"
struct nuc* get_nuc_by_num(struct rna* RNA, int nuc_num) {
  assert(RNA);
  if(RNA->first) {
    struct nuc* n = RNA->first;
    while(n) {
      if(nuc_num == n->num) {
	return n;
      }
      n = n->next;
    }
  }
  printf("No nuc found with number %d\n",nuc_num);
  exit(1);
}

// returns a pointer to an atom by searching for the serial number "atom_num"
struct atom* get_atom_by_num(struct nuc* n, int atom_num) {
  assert(n);
  if(n->first) {
    struct atom* a = n->first;
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

// retrieve the serial number of a nucletide 
int get_nuc_num(struct nuc* n) {
  assert(n);
  return n->num;
}

// retrieve the name of a nucleotide 
char* get_nuc_name(struct nuc* n) {
  assert(n);
  return n->name;
}

// retrieve the serial number of an atom
int get_atom_num(struct atom* a) {
  assert(a);
  return a->num;
}

// retrieve the name of an atom
char* get_atom_name(struct atom* a) {
  assert(a);
  return a->name;
}
