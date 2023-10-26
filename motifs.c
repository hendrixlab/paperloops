 /*
 * This library contains the definitions of structures and functions related to loop motifs
 * defining the lengths of each loop in a multiloop structure. e.g. [3,4,2,1]
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#include "motifs.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * data structures for the transform and for the library of transforms
 */

// this is just a dynamic array of "int"s, used to store primes and divisors
struct list {
  int* data;
  int len;
  int cap;
};

// this defines a single loop motif
struct motif {
  int* loop;
  int num;
  int id;
  struct motif* next;
};

// this is a group of loop motifs, used for recursively building libraries
struct motif_group {
  struct motif** data;
  int count; // current number of motifs
  int size; // should have a defined size
  int N; // number of loops in motif group
  int L; // total length of all loops in each motif
};

// this is the final output, a library of all motif types
struct motif_library {
  struct motif_group** data;
  int N_loops;
  int L_max;
};

// this is a "steptype" motif, translated from loop motif
// but is now a series of "bp" and "bb" transitions
struct steptypes {
  char** steps;
  int length;
};


// just a basic array of repeated "val" entries
int* create_array(int len, int val) {
  int* array = malloc(len * sizeof(int));
  for(int i=0; i<len; i++) {
    array[i] = val;
  }
  return array;
}

// this is creating a list, which is a dynamic array
struct list* create_list() {
  struct list* list = malloc(sizeof(list));
  list->cap = 1;
  list->len = 0;
  list->data = malloc((list->cap)*sizeof(int));
  return list;
}

// append to a list, update length as needed 
void list_append(struct list* list, int n) {
  if(list->len == list->cap) {
    int new_cap = 2*(list->cap);
    int* new_data = malloc(new_cap*sizeof(int));
    for(int i=0; i<list->len; i++) {
      new_data[i] = list->data[i];
    }
    list->data = new_data;
    list->cap = new_cap;
  }
  list->data[list->len] = n;
  list->len++;
}

// return the current length of a list/dynamic array
int list_length(struct list* list) {
  return list->len;
}

// print the contents of a list, for debugging
void print_list(struct list* list, char* title) {
  printf("%s: ",title);
  for(int i=0; i<list->len; i++) {
    printf("%d ", list->data[i]);
  }
  printf("\n");
}

int compare( const void* a, const void* b) {
  int int_a = * ( (int*) a );
  int int_b = * ( (int*) b );

  if(int_a > int_b) {
    return 1;
  } else if(int_a < int_b) {
    return -1;
  } else {
    return 0;
  }
}

/*************************************************************
 * Mathematical functions needed for enumerating loop motifs *
 *************************************************************/

struct list* get_divisors(int n) {
  // return the divisors of n
  int* divides = create_array(n+1,0);

  int total = 0;
  for(int i=1; i<=n; i++) {
    if(n % i == 0) {
      divides[i] = 1;
      total++;
    }
  }
  struct list* divisors = create_list();
  for(int i=1; i<=n; i++) {
    if(divides[i]) {
      list_append(divisors,i);
    }
  }  
  free(divides);
  return divisors;
}

struct list* get_primes(int n) {
  // compute an array of primes from 2 to n.
  int* prime = create_array(n+1,1);

  for(int i=2; i<=sqrt(n); i++) {
    if(prime[i]) {
      for(int j=i; i*j<=n; j++) {
	prime[i*j] = 0;
      }
    }
  }
  struct list* primes = create_list();
  for(int i=2; i<=n; i++) {
    if(prime[i]) {
      list_append(primes,i);
    }
  }
  free(prime);
  return primes;
}

int mobius(int n, struct list* primes) {
  if(n == 1) {
    return 1;
  } else {
    int total_p = 0;
    for(int i=0; i<list_length(primes); i++) {
      int p = primes->data[i];
      if(n % p == 0) {
	if(n % (p*p) == 0) {
	  return 0;
	} else {
	  total_p++;
	}
      }
    }
    if(total_p % 2 == 0) {
      // even number of prime factors
      return 1;
    } else {
      // odd number of prime factors
      return -1;
    }    
  }
}

int binomial(int n, int k) {
  // recursive computation of binomial coefficient
  // based on: https://www.geeksforgeeks.org/binomial-coefficient-dp-9/
  if (k == 0 || k == n) {
    return 1;
  }
  if (k > n) {
    return 0;
  }
  return binomial(n-1, k-1) + binomial(n-1, k);
}

/*****************************************************
 * Tools for theoretical enumeration of loop motifs *
 *****************************************************/

// this is the M(n,L) function, counts loops of n branches that sum to L
int aperiodic_necklace_counts_fixed_length(int n, int L, struct list* divisors, struct list* primes) {
  int sum = 0;
  for(int i=0; i<list_length(divisors); i++) {
    int d = divisors->data[i];
    if((L*d) % n == 0) {
      sum += mobius(n/d, primes) * binomial(L*d/n+d-1,d-1);
    }
  }
  return sum/n;  
}

// this is the N(n,L) function, Counts all loops of n branches, periodic and not, that sum to L
int** necklace_counts_fixed_length(int N_loops, int L_max) {
  int** M = malloc(sizeof(int*) * (N_loops+1));
  int** N = malloc(sizeof(int*) * (N_loops+1));
  for(int n=0; n<=N_loops; n++) {
    M[n] = malloc(sizeof(int) * (L_max+1));
    N[n] = malloc(sizeof(int) * (L_max+1));
  }
  
  struct list* primes = get_primes(N_loops); // all the primes we need =< N_loops
  for(int n=1; n<=N_loops; n++) {
    struct list* divisors = get_divisors(n); // divisors of n
    for(int L=0; L<=L_max; L++) {
      M[n][L] = aperiodic_necklace_counts_fixed_length(n,L,divisors,primes);
      N[n][L] = 0; // initializing N[n,L] updated below
    }
  }
  for(int n=1; n<=N_loops; n++) {
    struct list* divisors = get_divisors(n);
    for(int L=0; L<=L_max; L++) {
      for(int i=0; i<list_length(divisors); i++) {
	int d = divisors->data[i];
	if((L*d) % n == 0) {
	  N[n][L] += M[d][L*d/n];
	}
      }
    }
  }
  return N;
}

/*
 * Tools for building/adding/removing to a motif
 */

struct motif* create_motif(int num) {
  struct motif* motif = malloc(sizeof(struct motif));
  motif->num = num;
  motif->loop = malloc(sizeof(int)*num);
  motif->next = NULL;
  return motif;
}

void initialize_motif(struct motif* motif, int loop) {
  set_motif_loop(motif,0,loop);
  set_motif_id(motif,loop);
}

void set_motif_loop(struct motif* motif, int idx, int loop) {
  assert(motif);
  motif->loop[idx] = loop;
}

void set_motif_id(struct motif* motif, int id) {
  assert(motif);
  motif->id = id;
}

int get_motif_num(struct motif* motif) {
  return motif->num;
}

int get_motif_loop(struct motif* motif, int idx) {
  assert(motif);
  return motif->loop[idx];
}

int get_motif_id(struct motif* motif) {
  return motif->id;
}

int get_motif_sum(struct motif* motif) {
  int sum = 0;
  int num = get_motif_num(motif);
  for(int i=0; i<num; i++) {
    sum += get_motif_loop(motif,i);
  }
  return sum;
}

int print_motif(struct motif* motif, char c) {
  int n = get_motif_num(motif);
  printf("motif %d: (",motif->id);
  for(int i=0; i<n; i++) {
    printf("%d",motif->loop[i]);
    if(i==n-1) {
      printf(")%c",c);
    } else {
      printf(", ");
    }
  }  
}

int get_new_id(int id, int new_num, int new_loop, int L_max) {
  int new_id = id*(L_max+1) + new_loop; // next place in base-(L_max+1) number  
}

int compute_motif_id(struct motif* motif, int L_max) {
  int new_id = 0;
  int num = get_motif_num(motif);
  for(int n=0; n<num; n++) {
    new_id += get_motif_loop(motif,n)*pow(L_max+1,num-1-n);
  }
  return new_id;
}

struct motif* expand_motif(struct motif* motif, int new_loop, int L_max) {
  // adds to a motif a new loop length
  int num = get_motif_num(motif);
  int id = get_motif_id(motif);

  int new_num = num+1;
  struct motif* new_motif = create_motif(new_num);
  int new_id = get_new_id(id,new_num,new_loop,L_max);
  set_motif_id(new_motif,new_id);


  for(int i=0;i<num;i++) { // copy the old motif
    int loop = get_motif_loop(motif,i);
    set_motif_loop(new_motif,i,loop);
  }
  set_motif_loop(new_motif,num,new_loop); // add the new loop at i = num = new_num-1
  
  return new_motif;
}

struct motif* rotate_motif(struct motif* motif, int L_max) {
  // adds to a motif a new loop length
  int num = get_motif_num(motif);

  struct motif* rot_motif = create_motif(num);
    
  for(int i=1;i<num;i++) { // copy the old motif
    int loop = get_motif_loop(motif,i);
    set_motif_loop(rot_motif,i-1,loop);
  }
  int loop = get_motif_loop(motif,0);
  set_motif_loop(rot_motif,num-1,loop);
  
  int new_id = compute_motif_id(rot_motif,L_max);
  set_motif_id(rot_motif,new_id);

  // important: building in linked list of rotations attached to main motif
  //motif->next = rot_motif;
  
  return rot_motif;
}

struct motif* get_next_motif(struct motif* motif) {
  assert(motif);
  return motif->next;
}

struct motif* set_next_motif(struct motif* motif, struct motif* next) {
  assert(motif);
  return motif->next = next;
}

/*
 * Tools for building/adding/removing elements from motif_group
 */

struct motif_group* create_motif_group(int size) {
  struct motif_group* motif_group = malloc(sizeof(struct motif_group));  
  motif_group->count = 0;  
  motif_group->size = size;
  motif_group->data = malloc(sizeof(struct motif*) * size);
  return motif_group;
}

void free_motif_group(struct motif_group* motif_group) {
  assert(motif_group);
  free(motif_group->data);
  free(motif_group);
}

int get_motif_group_size(struct motif_group* motif_group) {
  assert(motif_group);
  return motif_group->size;
}

int get_motif_group_count(struct motif_group* motif_group) {
  assert(motif_group);
  return motif_group->count;
}

struct motif* get_motif(struct motif_group* motif_group, int idx) {
  assert(motif_group);
  if((idx < motif_group->size && idx >= 0)) {
    return motif_group->data[idx];
  }
  return NULL;
}

int motif_empty(struct motif_group* motif_group, int idx) {  
  if(idx < motif_group->size && idx >= 0) {
    if(motif_group->data[idx]) {
      return 0;
    }
    return 1;
  } else {
    printf("error: index out of range: %d not within [%d,%d]\n",idx,0,motif_group->size);
    exit(1);
  }
}

void add_motif(struct motif_group* motif_group, struct motif* motif) {
  assert(motif_group);
  assert(motif);

  int idx = get_motif_id(motif);
  assert(idx < motif_group->size && idx >= 0);
  
  motif_group->data[idx] = motif;
  motif_group->count++;
}

void append_motif(struct motif_group* motif_group, struct motif* motif) {
  assert(motif_group);
  assert(motif);

  int len = motif_group->count;
  if(len == motif_group->size) {
    // resize the motif_group
    int new_size = 2 * motif_group->size;
    struct motif** new_data = malloc(new_size * sizeof(struct motif*));
    for(int i=0; i<motif_group->count; i++) {
      new_data[i] = motif_group->data[i];
    }
    free(motif_group->data);
    motif_group->data = new_data;
    motif_group->size = new_size;
    //printf("resizing to %d\n",new_size);
  } 
  motif_group->data[len] = motif;
  motif_group->count++;
}

struct motif_group*** create_motif_lib(int N_loops, int L_max) {
  struct motif_group*** motif_lib = malloc((N_loops+1)*sizeof(struct motif_group**));
  for(int n=0;n<=N_loops;n++) {
    motif_lib[n] = malloc((L_max+1)*sizeof(struct motif_group*));
    for(int L=0; L<=L_max; L++) {
      motif_lib[n][L] = create_motif_group(1);
    }
  }
  return motif_lib;
}

struct motif_group** create_motif_list(int N_loops, int L_max) {
  struct motif_group** motif_list = malloc((N_loops+1)*sizeof(struct motif_group*));  
  motif_list[0] = create_motif_group(0); // sort of a place holder, has size=0
  motif_list[1] = create_motif_group(L_max+1); // hairpin loops, should be L_max+1 in size: 0,1,...,L_max
  
  for(int L=0;L<=L_max;L++) {
    struct motif* motif = create_motif(1);
    initialize_motif(motif,L); // one-loop motifs, setting l_0 = L;
    add_motif(motif_list[1],motif);
  }
  return motif_list;
}

struct motif_group*** build_motif_list(int N_loops, int L_max) {
  // build motif library
  struct motif_group*** motif_lib = create_motif_lib(N_loops,L_max);
  // build motif_group structure, initialize with singleton motifs:
  struct motif_group** motif_list = create_motif_list(N_loops,L_max);  

  struct list* prev_motif_ids = create_list();
  for(int L=0;L<=L_max;L++) {
    list_append(prev_motif_ids,L);
  }
  
  // recursively build motif_lib for different lengths    
  for(int n=2;n<=N_loops;n++) { // create multi-junction loops 
    int new_size = 1 + L_max * pow(L_max+1,n-1); // 1 + largest possible id for n loops that add up to L_max or less
    motif_list[n] = create_motif_group(new_size);
    struct list* motif_ids = create_list();
    for(int i=0;i<list_length(prev_motif_ids);i++) { // id of n-1 group
      int id = prev_motif_ids->data[i];
      struct motif* motif = get_motif(motif_list[n-1],id);
      int sum = get_motif_sum(motif);
      for(int L=0;L<=L_max-sum;L++) {
	int new_id = get_new_id(id,n,L,L_max);
	if(motif_empty(motif_list[n],new_id)) { // only add the new motif if it is not filled yet.
	  int L_tot = sum + L;
	  struct motif* new_motif = expand_motif(motif,L,L_max);
	  add_motif(motif_list[n],new_motif);
	  list_append(motif_ids,new_id); // append the id of each motif
	  append_motif(motif_lib[n][L_tot],new_motif); // only append the representative rotation
	  for(int j=0; j<n; j++) { // apart from repeats, motif with n loops has n rotations
	    struct motif* prev_motif = new_motif; // init. prev_motif to lexicographically lowest motif
	    new_motif = rotate_motif(new_motif,L_max); // first rotation
	    int new_id = get_motif_id(new_motif);
	    if(motif_empty(motif_list[n],new_id)) { // only add if it is not visited/stored yet
	      add_motif(motif_list[n],new_motif);
	      list_append(motif_ids,new_id); // append the id of each motif
	      set_next_motif(prev_motif,new_motif); // attach the rotation 
	    }
	  }
	}
      }
    }
    qsort(motif_ids->data,list_length(motif_ids),sizeof(int),compare);
    prev_motif_ids = motif_ids;
  }
  return motif_lib;
}

/* 
// original, slightly slower method for build_motif_list 
struct motif_group*** build_motif_list(int N_loops, int L_max) {
  // build motif library
  struct motif_group*** motif_lib = create_motif_lib(N_loops,L_max);
  // build motif_group structure, initialize with singleton motifs:
  struct motif_group** motif_list = create_motif_list(N_loops,L_max);  

  // recursively build motif_lib for different lengths    
  for(int n=2;n<=N_loops;n++) { // create multi-junction loops 
    int new_size = 1 + L_max * pow(L_max+1,n-1); // 1 + largest possible id for n loops that add up to L_max or less
    motif_list[n] = create_motif_group(new_size);
    int size = get_motif_group_size(motif_list[n-1]); // size of n-1
    for(int id=0;id<size;id++) { // id of n-1 group
      if(!motif_empty(motif_list[n-1],id)) { // checking if valid id for the n-1 group 
	struct motif* motif = get_motif(motif_list[n-1],id);
	int sum = get_motif_sum(motif);
	for(int L=0;L<=L_max-sum;L++) {
	  int new_id = get_new_id(id,n,L,L_max);
	  if(motif_empty(motif_list[n],new_id)) { // only add the new motif if it is not filled yet.
	    int L_tot = sum + L;
	    struct motif* new_motif = expand_motif(motif,L,L_max);
	    add_motif(motif_list[n],new_motif);	    
	    append_motif(motif_lib[n][L_tot],new_motif);
	    for(int j=0; j<n; j++) { // apart from repeats, motif with n loops has n rotations
	      struct motif* prev_motif = new_motif;
	      new_motif = rotate_motif(new_motif,L_max);
	      int new_id = get_motif_id(new_motif);
	      if(motif_empty(motif_list[n],new_id)) { // only add if it is not visited/stored yet
		add_motif(motif_list[n],new_motif);
		set_next_motif(prev_motif,new_motif);
	      }
 	    }
	  }
	}
      }
    }
  }
  return motif_lib;
}
 */

void print_motif_lib_counts(struct motif_group*** motif_lib, int N_loops, int L_max) {
  for(int n=0; n<N_loops; n++) {
    for(int L=0; L<=L_max; L++) {
      int count = motif_lib[n][L]->count;
      printf("%d %d %d\n",n,L,count);
    }
  }
}

/***********************
 * Tools for Steptypes *
 ***********************/

struct steptypes* create_steptypes(int length) {
  struct steptypes* steptypes = malloc(sizeof(struct steptypes));
  steptypes->steps = malloc(sizeof(char*)*length);
  steptypes->length = length;
  return steptypes;  
}

int steptypes_length(struct steptypes* steptypes) {
  assert(steptypes);
  return steptypes->length;
}

char* steptypes_get(struct steptypes* steptypes, int idx) {
  assert(steptypes);
  return steptypes->steps[idx];
}

void steptypes_set(struct steptypes* steptypes, int idx, char* step) {
  assert(steptypes);
  steptypes->steps[idx] = step;
}

struct steptypes* steptypes_splice(struct steptypes* steptypes, int start, int end) {
  struct steptypes* substeps = create_steptypes(end-start+1); // positions are inclusive
  for(int i=start; i<=end; i++) {
    char* step = steptypes_get(steptypes,i);
    steptypes_set(substeps,i,step);
  }
  return substeps;
}

struct steptypes* motif_to_steptypes(struct motif* motif) {
  int num = get_motif_num(motif);
  int sum = get_motif_sum(motif);
  struct steptypes* steptypes = create_steptypes(num+sum);
  int count = 0;
  steptypes_set(steptypes,count,"bp");
  count++;
  for(int i=0; i<num; i++) {    
    for(int j=0; j<motif->loop[i]; j++) {
      steptypes_set(steptypes,count,"bb");
      count++;
    }
    if(i < num-1) {
      steptypes_set(steptypes,count,"bb");
      count++;
    }
  }
  return steptypes;
}
