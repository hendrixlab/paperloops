/*
 * This library contains the definitions of structures and functions for implementing
 * the RNA data structure
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "RNA_molecule.h"

/*
 *  Data Structure Definitions 
 *  includes nodes for linked lists describing:
 *     -sequence of nucleotides (linked list of "nuc" nodes   
 *     -collection of atoms. There are three linked lists of atoms:
 *         -backbone atoms
 *         -sugar atoms
 *         -base atoms
 */

struct rna {
  struct nuc* first_nuc;
  struct nuc* last_nuc;
  int length;
  char* seq;
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
  double* bbta[6]; // backbone torsion angles
  double* sta[5];  // sugar torsion angles 
  struct nuc* next;
};

// atom node
struct atom {
  char* name;
  double* r;
  int num;
  char* type;
  struct atom* next;
};

/* 
 * Node creation subroutines 
 */

struct rna* create_rna() {
  struct rna* RNA = malloc(sizeof(struct rna));
  assert(RNA);
  RNA->first_nuc = NULL;
  RNA->last_nuc = NULL;
  RNA->length = 0;
  return RNA;
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
  return a;
}

/* 
 * Memory management/freeing subroutines 
 */

void free_rna(struct rna* RNA) {
  assert(RNA);
  
  while(!rna_is_empty(RNA)) {
    struct nuc* removed_nuc = RNA->first_nuc;
    RNA->first_nuc = removed_nuc->next;
    free_nuc(removed_nuc);
  }
  free(RNA);
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

int rna_is_empty(struct rna* RNA) {
  assert(RNA);
  return RNA->first_nuc == NULL;
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
 * Subroutines to add elements to the data structure 
 */

void add_nuc(struct rna* RNA, struct nuc* n) {
  if(RNA->first_nuc != NULL) { // if not the first nuc
    RNA->last_nuc->next = n;
    RNA->last_nuc = n;
  } else {  // the first (and currently last) nuc in the RNA chain
    RNA->first_nuc = n;
    RNA->last_nuc = n;
  }
  RNA->length++;
}

void add_atom(struct rna* RNA, struct atom* a) {
  assert(RNA);
  if(is_backbone_atom(a)) {
    // atom is part of backbone
    add_backbone_atom(RNA,a);
  } else if(is_sugar_atom(a)) {
    // atom is part of sugar
    add_sugar_atom(RNA,a);
  } else {
    //atom is part of the base
    add_base_atom(RNA,a);
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
  char* sugar_atoms[] = {"C2'","C1'","O4'"};
  int len = 3; // hard-coded 
  for(int i=0; i<len; i++) {
    if(strcmp(sugar_atoms[i],a->name) == 0) {
      return 1;
    }
  }
  return 0;
}

void add_backbone_atom(struct rna* RNA, struct atom* a) {
  a->type = "backbone";
  // first if there is a nucleotide
  if(RNA->first_nuc != NULL) {
    // next check if there is an atom attached to first nuc
    if(RNA->last_nuc->backbone_first != NULL) {
      RNA->last_nuc->backbone_last->next = a;
      RNA->last_nuc->backbone_last = a;
      RNA->last_nuc->size++;
    } else {
      // this is the first atom
      RNA->last_nuc->backbone_first = a;
      RNA->last_nuc->backbone_last = a;
    }
  } else {
    // error not expecting to add atom without nuc already in place
    printf("error: attempted to add an atom without a nuc in place.\n");
    exit(1);
  }
  RNA->last_nuc->size++;
}

void add_sugar_atom(struct rna* RNA, struct atom* a) {
  a->type = "sugar";
  // check if there is a nucleotide
  if(RNA->first_nuc != NULL) {
    // next check if there is an atom attached to first nuc
    if(RNA->last_nuc->sugar_first != NULL) {
      RNA->last_nuc->sugar_last->next = a;
      RNA->last_nuc->sugar_last = a;
      RNA->last_nuc->size++;
    } else {
      // this is the first atom
      RNA->last_nuc->sugar_first = a;
      RNA->last_nuc->sugar_last = a;
    }
  } else {
    // error not expecting to add atom without nuc already in place
    printf("error: attempted to add an atom without a nuc in place.\n");
    exit(1);
  }
  RNA->last_nuc->size++;
}

void add_base_atom(struct rna* RNA, struct atom* a) {
  a->type = "base";
  // check if there is a nucleotide
  if(RNA->first_nuc != NULL) {
    // next check if there is an atom attached to first nuc
    if(RNA->last_nuc->base_first != NULL) {
      RNA->last_nuc->base_last->next = a;
      RNA->last_nuc->base_last = a;
      RNA->last_nuc->size++;
    } else {
      // this is the first atom
      RNA->last_nuc->base_first = a;
      RNA->last_nuc->base_last = a;
    }
  } else {
    // error not expecting to add atom without nuc already in place
    printf("error: attempted to add an atom without a nuc in place.\n");
    exit(1);
  }
  RNA->last_nuc->size++;
}


/* 
 * Functions to access elements from RNA data structure 
 */

// this function checks if the last nucleotide in RNA has the number "nuc_num"
// used to avoid adding a nucleotide more than once 
int last_nuc_has_num(struct rna* RNA, int nuc_num) {
  assert(RNA);
  if(RNA->last_nuc != NULL) {
    return RNA->last_nuc->num == nuc_num;
  } else {
    return 0;
  }
}

// returns a pointer to the first nucleotide (nuc) in an RNA 
struct nuc* get_first_nuc(struct rna* RNA) {
  assert(RNA);
  if(RNA->first_nuc) {
    struct nuc* n = RNA->first_nuc;
    return n;
  }
  printf("Error: No nuc found for the RNA\n");
  exit(1);
}

// returns a pointer to a nuc by searching for the serial number "nuc_num"
struct nuc* get_nuc_by_num(struct rna* RNA, int nuc_num) {
  assert(RNA);
  if(RNA->first_nuc) {
    struct nuc* n = RNA->first_nuc;
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

//given a nuc return the next nuc in the RNA sequence
struct nuc* next_nuc(struct nuc* n) {
  assert(n);
  return n->next;
}

// returns a pointer to the previous backbone atom O3' relative to n
struct atom* get_prev_backbone_atom(struct nuc* n) {
  assert(n);
  if(n->backbone_first) {
    struct atom* a = n->backbone_first;
    return a;    
  } else {
    printf("Error: No backbone atom found.\n");
  }
}


// returns a pointer to the first backbone atom
struct atom* get_first_backbone_atom(struct nuc* n) {
  assert(n);
  if(n->backbone_first) {
    struct atom* a = n->backbone_first;
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
      a = a->next;
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

// retrieve the type of an atom
char* get_atom_type(struct atom* a) {
  assert(a);
  return a->type;
}

/* Functions for extracting the RNA sequence */

char* set_rna_sequence(struct rna* RNA) {
  assert(RNA);
  char* sequence; 
  if(RNA->length != 0) {
    sequence = malloc(sizeof(char) * RNA->length);
    struct nuc* n = RNA->first_nuc;
    strcpy(sequence,"");
    while(n != NULL) {
      strcat(sequence,n->name);
      n = n->next;
    }
  } else {
    printf("Error: RNA has length zero.\n");
    exit(0);
  }
  RNA->seq = sequence;
}

char* get_rna_sequence(struct rna* RNA) {
  assert(RNA);
  if(RNA->seq == NULL) {
    set_rna_sequence(RNA);
  }
  return RNA->seq;
}
