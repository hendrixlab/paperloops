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
#include "vectorfun.h"
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
  struct nuc* sentinel;
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
};

/* 
 * Node creation subroutines 
 */

struct rna* create_rna() {
  struct rna* RNA = malloc(sizeof(struct rna));
  RNA->sentinel = malloc(sizeof(struct nuc));
  RNA->sentinel->num = 0;
  RNA->sentinel->name = "sentinel";
  RNA->sentinel->next = RNA->sentinel;
  RNA->sentinel->prev = RNA->sentinel;
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
  return a;
}

/* 
 * Memory management/freeing subroutines 
 */

void free_rna(struct rna* RNA) {
  assert(RNA);
  
  while(!rna_is_empty(RNA)) {
    struct nuc* removed_nuc = get_first_nuc(RNA);
    RNA->sentinel = removed_nuc->next;
    removed_nuc->next->prev = RNA->sentinel;
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
  return RNA->sentinel->next == RNA->sentinel;
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
  RNA->sentinel->prev->next = n;
  n->prev = RNA->sentinel->prev;
  RNA->sentinel->prev = n;
  n->next = RNA->sentinel;
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
  if(rna_is_empty(RNA)) {   // first if there is a nucleotide
    printf("error: attempted to add a backbone atom without a nuc in place.\n");
    exit(1);
  } else {  // next check if there is an atom attached to first nuc
    if(RNA->sentinel->prev->backbone_first != NULL) {
      RNA->sentinel->prev->backbone_last->next = a;
      a->prev = RNA->sentinel->prev->backbone_last;
      RNA->sentinel->prev->backbone_last = a;
      RNA->sentinel->prev->size++;
    } else {   // this is the first atom
      RNA->sentinel->prev->backbone_first = a;
      RNA->sentinel->prev->backbone_last = a;
    }
  }
  RNA->sentinel->prev->size++;
}

void add_sugar_atom(struct rna* RNA, struct atom* a) {
  a->type = "sugar";
  if(rna_is_empty(RNA)) {  // check if there is a nucleotide
    printf("error: attempted to add a sugar atom without a nuc in place.\n");
    exit(1);
  } else {  // next check if there is an atom attached to first nuc
    if(RNA->sentinel->prev->sugar_first != NULL) {
      RNA->sentinel->prev->sugar_last->next = a;
      a->prev = RNA->sentinel->prev->sugar_last;
      RNA->sentinel->prev->sugar_last = a;
      RNA->sentinel->prev->size++;
    } else {
      // this is the first atom
      RNA->sentinel->prev->sugar_first = a;
      RNA->sentinel->prev->sugar_last = a;
    }
  }
  RNA->sentinel->prev->size++;
}

void add_base_atom(struct rna* RNA, struct atom* a) {
  a->type = "base";
  if(rna_is_empty(RNA)) { // check if there is a nucleotide
    printf("error: attempted to add a base atom without a nuc in place.\n");
    exit(1);
  } else {
    // next check if there is an atom attached to first nuc
    if(RNA->sentinel->prev->base_first != NULL) {
      RNA->sentinel->prev->base_last->next = a;
      a->prev = RNA->sentinel->prev->base_last;
      RNA->sentinel->prev->base_last = a;
      RNA->sentinel->prev->size++;
    } else {
      // this is the first atom
      RNA->sentinel->prev->base_first = a;
      RNA->sentinel->prev->base_last = a;
    }
  }
  RNA->sentinel->prev->size++;
}

/******************************************************** 
 * Functions to access elements from RNA data structure *
 ********************************************************/

// this function checks if the last nucleotide in RNA has the number "nuc_num"
// used to avoid adding a nucleotide more than once 
int new_nuc_num(struct rna* RNA, int nuc_num) {
  assert(RNA);
  if(rna_is_empty(RNA)) {
    return 1;
  } else {
    return RNA->sentinel->prev->num != nuc_num;
  }
}

// returns a pointer to the first nucleotide (nuc) in an RNA 
struct nuc* get_first_nuc(struct rna* RNA) {
  assert(RNA);
  if(!rna_is_empty(RNA)) {
    struct nuc* n = RNA->sentinel->next;
    return n;
  }
  printf("Error: No nuc found for the RNA\n");
  exit(1);
}

// returns a pointer to a nuc by searching for the serial number "nuc_num"
struct nuc* get_nuc_by_num(struct rna* RNA, int nuc_num) {
  assert(RNA);
  if(rna_is_empty(RNA)) {
    struct nuc* n = RNA->sentinel->next;
    while(before_end(RNA,n)) {
      if(nuc_num == n->num) {
	return n;
      }
      n = next_nuc(RNA,n);
    }
  }
  printf("No nuc found with number %d\n",nuc_num);
  exit(1);
}

// checks if the nuc "n" is within the RNA (not the sentinel)
int before_end(struct rna* RNA, struct nuc* n) {
  assert(RNA);
  assert(n);
  return n != RNA->sentinel;
}

//given a nuc return the next nuc in the RNA sequence
struct nuc* next_nuc(struct rna* RNA, struct nuc* n) {
  assert(n);
  return n->next;
}

//given a nuc return the next nuc in the RNA sequence
struct nuc* prev_nuc(struct rna* RNA, struct nuc* n) {
  assert(n);
  return n->prev;
}

// returns a pointer to the previous backbone atom O3' relative to n
struct atom* prev_nuc_last_backbone_atom(struct rna* RNA, struct nuc* n) {
  assert(RNA);
  assert(n);
  if(!before_end(RNA,prev_nuc(RNA,n))) {
    return NULL;
  } else {
    return n->prev->backbone_last;
  }
}

// returns a pointer to the next backbone atom P relative to n+1
struct atom* next_nuc_first_backbone_atom(struct rna* RNA, struct nuc* n) {
  assert(RNA);
  assert(n);
  if(!before_end(RNA,next_nuc(RNA,n))) {
    return NULL;
  } else {
    return n->next->backbone_first;
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

// returns a pointer to the next atom in a set given pointer to atom
// e.g. the next backbone atom
struct atom* prev_atom(struct atom* a) {
  assert(a);
  return a->prev;
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

/******************************
 * Torsion angle calculations *
 ******************************/

double torsion(struct atom* p1, struct atom* p2, struct atom* p3, struct atom *p4) {
  printf("p1 - %f,%f,%f\n",p1->r[0],p1->r[1],p1->r[2]);
  printf("p2 - %f,%f,%f\n",p2->r[0],p2->r[1],p2->r[2]);
  printf("p3 - %f,%f,%f\n",p3->r[0],p3->r[1],p3->r[2]);
  printf("p4 - %f,%f,%f\n",p4->r[0],p4->r[1],p4->r[2]);
  double* vec_a = malloc(3 * sizeof(double));
  double* vec_b = malloc(3 * sizeof(double));
  double* vec_c = malloc(3 * sizeof(double));
  vec_a = subtract(p2->r,p1->r);
  vec_b = subtract(p3->r,p2->r);
  vec_c = subtract(p4->r,p3->r);
  printf("vec_a - %f,%f,%f\n",vec_a[0],vec_a[1],vec_a[2]);
  printf("vec_b - %f,%f,%f\n",vec_b[0],vec_b[1],vec_b[2]);
  printf("vec_c - %f,%f,%f\n",vec_c[0],vec_c[1],vec_c[2]);
  //double dot_prod = (dot(vec_b,vec_c) * dot(vec_a,vec_b))/dot(vec_b,vec_b) - dot(vec_a,vec_c);  
  double* proj_a = orthonormal_projection(vec_a,vec_b);
  double* proj_c = orthonormal_projection(vec_c,vec_a);
  printf("proj_a - %f,%f,%f\n",proj_a[0],proj_a[1],proj_a[2]);
  printf("proj_c - %f,%f,%f\n",proj_c[0],proj_c[1],proj_c[2]);
  double angle = get_angle(proj_a,proj_c);
  printf("angle: %f\n",angle);
  double sign = dot(cross(proj_a,proj_c),vec_b);
  if(sign < 0) {
    angle = -angle;
  }
  return angle;
}

void compute_torsion_angles(struct rna* RNA) {
  struct nuc* n = get_first_nuc(RNA);
  int FIRST = 1;
  while(before_end(RNA,n)) {
    struct atom* p1 = prev_nuc_last_backbone_atom(RNA,n);
    struct atom* p2 = get_first_backbone_atom(n);
    struct atom* p3 = next_atom(p2);
    struct atom* p4 = next_atom(p3);
    int count;
    if(strcmp(get_atom_name(p2),"P") == 0) {
      count = 0;
    } else if(strcmp(get_atom_name(p2),"O5'") == 0) {
      count = 1;
    } else {
      printf("Unexpected first atom in backbone: \"%s\"\n",get_atom_name(p2));
      exit(1);
    }
    int CONTINUE = 1;
    int NEXT_NUC = 0;
    while(CONTINUE) {
      double angle = -100;
      if(p1 && p2 && p3 && p4) {
	angle = torsion(p1,p2,p3,p4);
      }
      printf("%d - %s: ",get_nuc_num(n),get_nuc_name(n));
      char* label1 = p1 ? p1->name : "None";
      char* label2 = p2 ? p2->name : "None";
      char* label3 = p3 ? p3->name : "None";
      char* label4 = p4 ? p4->name : "None";
      printf("%d - %f (%s,%s,%s,%s) %d\n",count,angle,label1,label2,label3,label4,NEXT_NUC);
      // update to the next backbone atom quadruple
      p1 = p2;
      p2 = p3;
      p3 = p4;
      if(p4 && next_atom(p4)) {
	p4 = next_atom(p4);
	//printf("retrieved atom method 1\n");
	if(NEXT_NUC > 1) {
	  CONTINUE = 0; // allow two atoms from next nucleotide P & O5'
	} else if(NEXT_NUC > 0) {
	  NEXT_NUC++;
	}
      } else {
	if(before_end(RNA,next_nuc(RNA,n))) {
	  p4 = next_nuc_first_backbone_atom(RNA,n);
	  NEXT_NUC++;
	  //printf("retrieved atom method 2\n");
	} else {
	  CONTINUE = 0;
	}
      }
      count++;
    }
    printf("\n");
    n = next_nuc(RNA,n);
  }
}

/* Functions for extracting the RNA sequence */

char* set_rna_sequence(struct rna* RNA) {
  assert(RNA);
  char* sequence; 
  if(RNA->length != 0) {
    sequence = malloc(sizeof(char) * RNA->length);
    struct nuc* n = RNA->sentinel->next;
    strcpy(sequence,"");
    while(before_end(RNA,n)) {
      strcat(sequence,n->name);
      n = next_nuc(RNA,n);
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
  } else if(strlen(RNA->seq) != RNA->length) {
    set_rna_sequence(RNA);    
  }
  return RNA->seq;
}
