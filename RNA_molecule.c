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
  a->sugar_next = NULL;
  a->sugar_prev = NULL;
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
  //printf("adding: %s\n",a->name);
  if(is_backbone_atom(a)) {
    // atom is part of backbone
    //printf("adding to backbone\n");
    add_backbone_atom(RNA,a);
  }
  if(is_sugar_atom(a)) {
    // atom is part of sugar
    //printf("adding to sugar\n");
    add_sugar_atom(RNA,a);
  }
  if((!is_backbone_atom(a))&&(!is_sugar_atom(a))) {
    //atom is part of the base
    //printf("adding to base\n");
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
  char* sugar_atoms[] = {"C1'","O4'","C4'","C3'","C2'"};
  int len = 5; // hard-coded 
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
    // **note: sugar atoms use separate pointers 
    if(RNA->sentinel->prev->sugar_first != NULL) {
      RNA->sentinel->prev->sugar_last->sugar_next = a;
      a->sugar_prev = RNA->sentinel->prev->sugar_last;
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
    return RNA->sentinel->prev->num != nuc_num; // return true if nuc_num is a new nuc_num
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
struct atom* print_sugar_atoms(struct nuc* n) {
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
  printf("No atom found with name: %s\n",atom_name);
  exit(1);
}

void sort_sugar_atoms(struct nuc* n) {
  struct atom* new_sugar_first = NULL;
  struct atom* new_sugar_last = NULL;
  char* sugar_atoms[] = {"C1'","O4'","C4'","C3'","C2'"};
  int len = 5; // hard-coded
  for(int i=0; i<len; i++) {
    //print_sugar_atoms(n);
    //printf("removing %d %s\n",i,sugar_atoms[i]);
    struct atom* a = remove_sugar_atom_by_name(n,sugar_atoms[i]);
    //printf("sort: adding in %s\n",a->name);
    //print_sugar_atoms(n);
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

double torsion1(struct atom* p1, struct atom* p2, struct atom* p3, struct atom *p4) {
  // method 1: most intutive method, from: https://x3dna.org/highlights/how-to-calculate-torsion-angle
  if(p1 && p2 && p3 && p4) {
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

void compute_backbone_torsion(struct rna* RNA) {
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
      n->bbta[0] = NAN;
      count = 1;
    } else {
      printf("Unexpected first atom in backbone: \"%s\"\n",get_atom_name(p2));
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
	if(before_end(RNA,next_nuc(RNA,n))) {
	  p4 = next_nuc_first_backbone_atom(RNA,n);
	  NEXT_NUC++;
	} else {
	  CONTINUE = 0;
	}
      }
      count++;
    }
    n = next_nuc(RNA,n);
  }
}

void compute_sugar_torsion(struct rna* RNA) {
  struct nuc* n = get_first_nuc(RNA);  
  while(before_end(RNA,n)) {
    // sort atoms in bonded order
    sort_sugar_atoms(n);
    int count = 0;
    struct atom* p2 = get_first_sugar_atom(n);
    struct atom* p1 = prev_sugar_atom(p2);
    struct atom* p3 = next_sugar_atom(p2);
    struct atom* p4 = next_sugar_atom(p3);;
    while(count<5) {
      n->sta[count] = torsion1(p1,p2,p3,p4); // store torsion angle    
      printf("%d %s: ",n->num,n->name);
      printf("computing torsion for: %s %s %s %s -> ",p1->name,p2->name,p3->name,p4->name);
      printf("%.2f\n",n->sta[count]);
      p1 = p2;
      p2 = p3;
      p3 = p4;
      p4 = next_sugar_atom(p3);
      count++;
    }
    n = next_nuc(RNA,n);
  }
}


void print_backbone_torsion(struct rna* RNA) {
  char* angle_name[] = {"alpha","beta","gamma","delta","epsilon","zeta"};
  struct nuc* n = get_first_nuc(RNA);
  int FIRST = 1;
  /* header */
  printf("pos\tnt");
  for(int i=0; i<6; i++) {
    printf("\t%s",angle_name[i]);
  }
  printf("\n");
  /* data table */
  while(before_end(RNA,n)) {
    printf("%d\t%s",n->num,n->name);
    for(int i=0; i<6; i++) {
      printf("\t%.2f",n->bbta[i]);
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
