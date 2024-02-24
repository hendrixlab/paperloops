/*
 * This library contains the definitions of structures and functions related to loop motifs
 * defining the lengths of each loop in a multiloop structure. e.g. [3,4,2,1]
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 * In collaboration with Rhiju Das, Stanford
 * translated/derived from toyfold3 with permission: https://github.com/rhiju/toyfold3 (C) Rhiju Das
 */

#include "motifs.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * data structures for the motifs and combinatorics
 */

// this is just a dynamic array of "int"s, used to store primes and divisors
struct list {
  long int* data;
  int len;
  int cap;
};

// node for linked list to avoid hash table collisions
struct hash_node {
  long int key;
  void* value;
  struct hash_node* next;
};

// a hash table structure
struct hash_table {
  struct hash_node** array;
  int size;
  int total;
};

// this defines a single loop motif
struct motif {
  int* loop;
  int num;
  long int id;
  struct motif* next;
};

// this is a group of loop motifs, used for recursively building libraries
struct motif_group {
  struct motif** data; // array of pointers to motifs
  int count; // current number of motifs
  int size; // should have a defined size/capacity >= count
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
// but is now a series of "BP" and "BB" transitions
struct steptypes {
  char** steps;
  int length;
};

// just a basic array of repeated "val" entries
int* create_rep_array(int len, int val) {
  int* array = malloc(len * sizeof(int));
  for(int i=0; i<len; i++) {
    array[i] = val;
  }
  return array;
}

// this is creating a list, which is a dynamic array
struct list* create_list() {
  struct list* L = malloc( sizeof(struct list) );
  L->cap = 1;
  L->len = 0;
  L->data = malloc((L->cap) * sizeof(long int));
  return L;
}

void free_list(struct list* list) {
  free(list->data);
  free(list);
}

// append to a list, update length as needed 
void list_append(struct list* list, long int n) {
  if(list->len == list->cap) {
    int new_cap = 2*(list->cap);
    long int* new_data = malloc(new_cap*sizeof(long int));
    for(int i=0; i<list->len; i++) {
      new_data[i] = list->data[i];
    }
    free(list->data);
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
    printf("%l ", list->data[i]);
  }
  printf("\n");
}

/************************
 * hash table functions *
 ************************/

int hash_function(struct hash_table* hash_table, long int key) {
  return key % hash_table->size;
}

struct hash_table* hash_table_create(int array_size) {
  struct hash_table* hash_table = malloc(sizeof(struct hash_table));
  assert(hash_table);
  hash_table->total = 0;
  hash_table->size = array_size;

  // allocate memory for the array, assign lists to NULL 
  hash_table->array = malloc(array_size * sizeof(struct hash_node*));

  // assign each term to NULL
  for(int i=0; i < hash_table->size; i++) {
    hash_table->array[i] = NULL;
  }

  // return the hash_table
  return hash_table;
}

void hash_table_free(struct hash_table* hash_table) {
  assert(hash_table);
  /*
   * loop through the array, and delete each node
   */
  for(int i=0; i < hash_table->size; i++) {  
    struct hash_node* current = hash_table->array[i];
    while (current != NULL) {
      hash_table->array[i] = current->next;
      assert(current);
      free(current);
      current  = hash_table->array[i];
    }
  }
  free(hash_table);
}

void hash_table_reset(struct hash_table* hash_table) {
  assert(hash_table);
  /*
   * loop through the array, and delete each node
   */
  for(int i=0; i < hash_table->size; i++) {  
    struct hash_node* current = hash_table->array[i];
    while (current != NULL) {
      hash_table->array[i] = current->next;
      assert(current);
      free(current);
      current = hash_table->array[i];
      // reduce the total by 1 for each removed node
      hash_table->total--;
    }
  }
  // at the point, the hash_table should be an array with NULL pointers only
}

void hash_table_set(struct hash_table* hash_table, long int key, void* value) {
  assert(hash_table);
  /*
   * Fill out the new node at put it at the beginning of the linked list 
   * for the hash index, defined using the key/name of the product 
   */
  struct hash_node* new_node = malloc(sizeof(struct hash_node));
  assert(new_node);

  new_node->key = key;
  new_node->value = value;
  new_node->next = NULL;

  int hash_index = hash_function(hash_table, key);  

  // add the node to the beginning of the list 
  new_node->next = hash_table->array[hash_index];
  hash_table->array[hash_index] = new_node;

  // update the total to reflect new addition
  hash_table->total++;
}

void* hash_table_get(struct hash_table* hash_table, long int key) {
  assert(hash_table);
  int hash_index = hash_function(hash_table, key);  

  // add the node to the beginning of the list 
  struct hash_node* node = hash_table->array[hash_index];
  while(node) {
    if(key == node->key) {
      return node->value;
    }
    node = node->next;
  }
  return NULL;
}

int hash_table_empty(struct hash_table* hash_table, long int key) {
  assert(hash_table);
  int hash_index = hash_function(hash_table, key);  

  // add the node to the beginning of the list 
  struct hash_node* node = hash_table->array[hash_index];
  while(node) {
    if(key == node->key) {
      // key not found, return False
      return 0;
    }
    node = node->next;
  }
  // key not found, return True
  return 1;
}

float hash_table_load(struct hash_table* hash_table) {
  int num_col = 0;
  
  // loop through the array, then traverse the linked list
  for(int i=0; i < hash_table->size; i++) {  
    if(hash_table->array[i]) {
      struct hash_node* current = hash_table->array[i];
      while (current->next) {
        current = current->next;
        num_col++;
      }
    }
  }
  return (float)num_col / (float)hash_table->size;
}

/****************************
 * end hash table functions *
 ****************************/

int compare( const void* a, const void* b) {
  long int int_a = * ( (long int*) a );
  long int int_b = * ( (long int*) b );

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
  int* divides = create_rep_array(n+1,0);
  
  for(int i=1; i<=n; i++) {
    if(n % i == 0) {
      divides[i] = 1;
    }
  }
  struct list* divisors = create_list();
  for(int i=1; i<=n; i++) {
    if(divides[i]) {
      list_append(divisors,(long int)i);
    }
  }  
  free(divides);
  return divisors;
}

struct list* get_primes(int n) {
  // compute an array of primes from 2 to n.
  int* prime = create_rep_array(n+1,1);

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
      long int p = primes->data[i];
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

int aperiodic_necklace_counts_fixed_length(int n, int L, struct list* divisors, struct list* primes) {
  // this is the M(n,L) function, counts loops of n branches that sum to L that are equivalent under rotation
  int sum = 0;
  for(int i=0; i<list_length(divisors); i++) {
    long int d = divisors->data[i];
    if((L*d) % n == 0) {
      sum += mobius(n/d, primes) * binomial(L*d/n+d-1,d-1);
    }
  }
  return sum/n;  
}

int** create_necklace_counts(int N_loops, int L_max) {
  int** N = malloc(sizeof(int*) * (N_loops+1));
  for(int n=0; n<=N_loops; n++) {
    N[n] = malloc(sizeof(int) * (L_max+1));
  }
  return N;
}

void free_necklace_counts(int** N, int N_loops) {
  assert(N);
  for(int n=0; n<=N_loops; n++) {
    free(N[n]);
  }
  free(N);  
}

int** necklace_counts_fixed_length(int N_loops, int L_max) {
  // this is the N(n,L) function, counts all loops of n branches, periodic and not, that sum to L
  int** N = create_necklace_counts(N_loops, L_max);
  int** M = create_necklace_counts(N_loops, L_max);
  
  struct list* primes = get_primes(N_loops); // all the primes we need =< N_loops
  for(int n=1; n<=N_loops; n++) {
    struct list* divisors = get_divisors(n); // divisors of n
    for(int L=0; L<=L_max; L++) {
      M[n][L] = aperiodic_necklace_counts_fixed_length(n,L,divisors,primes);
      N[n][L] = 0; // initializing N[n,L] updated below
    }
    free_list(divisors);
  }
  for(int n=1; n<=N_loops; n++) {
    struct list* divisors = get_divisors(n);
    for(int L=0; L<=L_max; L++) {
      for(int i=0; i<list_length(divisors); i++) {
	long int d = divisors->data[i];
	if((L*d) % n == 0) {
	  // all necklace counts are the sum of valid periodic and aperiodic counts 
	  N[n][L] += M[d][L*d/n];
	}
      }
    }
    free_list(divisors);
  }
  // clean up memory
  free_list(primes);
  free_necklace_counts(M,N_loops);
  // return N where N[n][L] = number of necklaces with n loops totalling L
  return N;
}

/*************************************************
 * Tools for building/adding/removing to a motif *
 *************************************************/

struct motif* create_motif(int num) {
  struct motif* motif = malloc(sizeof(struct motif));
  motif->num = num;
  motif->loop = malloc(sizeof(int)*num);
  motif->next = NULL;
  return motif;
}

void free_motif(struct motif* motif) {
  assert(motif);
  free(motif->loop);
  free(motif);
}

void initialize_motif(struct motif* motif, int loop) {
  // initialize a singleton n=1 motif. only works for these
  set_motif_loop(motif,0,loop);
  set_motif_id(motif,loop);
}

void set_motif_loop(struct motif* motif, int idx, int loop) {
  assert(motif);
  motif->loop[idx] = loop;
}

void set_motif_id(struct motif* motif, long int id) {
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

long int get_motif_id(struct motif* motif) {
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

void print_motif(struct motif* motif, char c) {
  int n = get_motif_num(motif);
  //printf("%d\t(",motif->id);
  for(int i=0; i<n; i++) {
    printf("%d",motif->loop[i]);
    if(i==n-1) {
      printf("%c",c);
    } else {
      printf("-");
    }
  }  
}

void fprint_motif(FILE* F, struct motif* motif, char c) {
  int n = get_motif_num(motif);
  //fprintf(F,"%d\t(",motif->id);
  for(int i=0; i<n; i++) {
    fprintf(F,"%d",motif->loop[i]);
    if(i==n-1) {
      fprintf(F,"%c",c);
    } else {
      fprintf(F,"-");
    }
  }  
}

long int get_new_id(long int id, int new_num, int new_loop, int L_max) {
  long int new_id = id*(L_max+1) + new_loop; // next place in base-(L_max+1) number
  return new_id;
}

long int compute_motif_id(struct motif* motif, int L_max) {
  long int new_id = 0;
  int num = get_motif_num(motif);
  for(int n=0; n<num; n++) {
    new_id += get_motif_loop(motif,n)*pow(L_max+1,num-1-n);
  }
  return new_id;
}

struct motif* expand_motif(struct motif* motif, int new_loop, int L_max) {
  // adds to a motif a new loop length
  int num = get_motif_num(motif);
  long int id = get_motif_id(motif);

  int new_num = num+1;
  struct motif* new_motif = create_motif(new_num);
  long int new_id = get_new_id(id,new_num,new_loop,L_max);
  set_motif_id(new_motif,new_id);

  for(int i=0;i<num;i++) { // copy the old motif
    int loop = get_motif_loop(motif,i);
    set_motif_loop(new_motif,i,loop);
  }
  set_motif_loop(new_motif,num,new_loop); // add the new loop at i = num = new_num-1
  
  return new_motif;
}

struct motif* rotate_motif(struct motif* motif, int L_max) {
  // rotates motif, x0-x1-x2-...-xn -> x1-x2-...-xn-x0
  // can not free incoming motif here because stored in build_motif_lib(). Must free at end of program run
  int num = get_motif_num(motif);

  struct motif* rot_motif = create_motif(num);
    
  for(int i=1;i<num;i++) { // copy the old motif
    int loop = get_motif_loop(motif,i);
    set_motif_loop(rot_motif,i-1,loop);
  }
  int loop = get_motif_loop(motif,0);
  set_motif_loop(rot_motif,num-1,loop);
  
  long int new_id = compute_motif_id(rot_motif,L_max);
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
  motif_group->data = calloc(size, sizeof(struct motif*)); // array of pointers to motifs
  //for(int i=0;i<size;i++) motif_group->data[i] = NULL; 
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
  if(motif_group) {
    return motif_group->count;
  } else {
    return 0;
  }
}

struct motif* get_motif(struct motif_group* motif_group, int idx) {
  // fix: needs to be hash_table for all applications? 
  assert(motif_group);
  if((idx < motif_group->size && idx >= 0)) {
    return motif_group->data[idx];
  }
  return NULL;
}

int motif_empty(struct motif_group* motif_group, int idx) {  
  // might be defunct
  if(0 <= idx && idx < motif_group->size) {
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
  motif_group->data[len] = motif; // add to last spot, this need not be motif->id
  motif_group->count++;
}

struct motif_group*** create_motif_lib(int N_loops, int L_max) {
  // motif_group array keyed by both number of loops n and the total length L
  struct motif_group*** motif_lib = malloc((N_loops+1)*sizeof(struct motif_group**));
  for(int n=0; n<=N_loops; n++) {
    motif_lib[n] = malloc((L_max+1)*sizeof(struct motif_group*));
    for(int L=0; L<=L_max; L++) {
      motif_lib[n][L] = create_motif_group(1);
    }
  }
  return motif_lib;
}

struct hash_table** create_motif_list(int N_loops, int L_max) {
  struct hash_table** motif_list = malloc((N_loops+1)*sizeof(struct hash_table*));  
  motif_list[0] = hash_table_create(0); // sort of a place holder, has size=0
  motif_list[1] = hash_table_create(L_max+1); // hairpin loops, should be L_max+1 in size: 0,1,...,L_max
  
  return motif_list;
}

void free_motif_list(struct hash_table** motif_list, int N_loops) {
  // motif list is an array of motif groups, indexed by the number of loops "n"
  for(int n=0; n<=N_loops; n++) {
    hash_table_free(motif_list[n]);
  }
  free(motif_list);
}

void free_motif_lib(struct motif_group*** motif_lib, int N_loops, int L_max) {
  // a motif lib is a 2D array of motif_groups, indexed by n and L
  assert(motif_lib);
  for(int n=0; n<=N_loops; n++) { // number of loops n
    for(int L=0; L<=L_max; L++) { // total length of loops L
      // here is the only place we explicitly go through and free the individual motifs
      int count = get_motif_group_count(motif_lib[n][L]);
      for(int i=0; i<count; i++) {
	struct motif* motif = get_motif(motif_lib[n][L],i);
	while(motif) {
	  struct motif* next_motif = get_next_motif(motif);
	  free_motif(motif);
	  motif = next_motif;
	}
	//free_motif(motif);

	/*
	if(motif) {
	  while(get_next_motif(motif)) {
	    struct motif* next_motif = get_next_motif(motif);
	    free_motif(motif);
	    motif = next_motif;
	  }
	  free_motif(motif);
	}
	*/
      }
      free_motif_group(motif_lib[n][L]);
    }
    free(motif_lib[n]);
  }
  free(motif_lib);
}

struct motif_group*** build_motif_lib(int N_loops, int L_max, int** N) {
  // create data structures
  struct hash_table** motif_list = create_motif_list(N_loops,L_max); // indexed by n only, all motifs for that n
  struct motif_group*** motif_lib = create_motif_lib(N_loops,L_max); // indexed by n and L, only "lowest" rotation
  struct list* prev_motif_ids = create_list();

  // initialize with singleton motifs:
  for(int L=0;L<=L_max;L++) {
    struct motif* motif = create_motif(1);
    initialize_motif(motif,L); // one-loop motifs, setting l_0 = L;
    hash_table_set(motif_list[1],L,motif);
    append_motif(motif_lib[1][L],motif);
    list_append(prev_motif_ids,L);
  }
  
  // recursively build motif_lib for different lengths    
  for(int n=2;n<=N_loops;n++) { // create multi-junction loops 
    //int new_size = 1 + L_max * pow(L_max+1,n-1); // 1 + largest possible id for n loops that add up to L_max or less
    int new_size = N[n][L_max];
    printf("New motif size for n=%d: %d\n",n,new_size);
    motif_list[n] = hash_table_create(new_size);
    struct list* motif_ids = create_list();
    for(int i=0;i<list_length(prev_motif_ids);i++) { // id of n-1 group
      long int id = prev_motif_ids->data[i];
      // get the motif, cast to a struct motif* pointer
      struct motif* motif = (struct motif*) hash_table_get(motif_list[n-1],id);
      int sum = get_motif_sum(motif);
      for(int L=0;L<=L_max-sum;L++) {
	long int new_id = get_new_id(id,n,L,L_max);
	if(new_id < 0) {
	  printf("%l\t%d\t%d\n",new_id,n,L);
	}
	if(hash_table_empty(motif_list[n],new_id)) { // only add the new motif if it is not filled yet.
	  int L_tot = sum + L;
	  struct motif* new_motif = expand_motif(motif,L,L_max);
	  hash_table_set(motif_list[n],new_id,new_motif);
	  append_motif(motif_lib[n][L_tot],new_motif); //lexicographically lowest, only one added to lib
 	  list_append(motif_ids,new_id); // append the id of each motif
	  for(int j=0; j<n; j++) { // apart from repeats, motif with n loops has n rotations
	    struct motif* prev_motif = new_motif; // init. prev_motif to lexicographically lowest motif
	    new_motif = rotate_motif(new_motif,L_max); // first rotation
	    long int new_id = get_motif_id(new_motif);
	    if(new_id < 0) {
	      printf("%l\t%d\t%d\n",new_id,n,L);
	    }
	    if(hash_table_empty(motif_list[n],new_id)) { // only add if it is not visited/stored yet
	      hash_table_set(motif_list[n],new_id,new_motif); // each rotation gets added to motif_list
	      list_append(motif_ids,new_id); // append the id of each rotation to motif_ids
	      set_next_motif(prev_motif,new_motif); // attach the rotation so motif_lib has all of them
	    }
	  }
	}
      }
    }
    qsort(motif_ids->data,list_length(motif_ids),sizeof(long int),compare);
    free_list(prev_motif_ids); // this should free memory allocated for motif_ids b/c of next line
    prev_motif_ids = motif_ids;
  }
  // final memory clean-up
  free_list(prev_motif_ids);
  free_motif_list(motif_list,N_loops);
  return motif_lib;
}


/*
struct motif_group*** build_motif_lib(int N_loops, int L_max) {
  // create data structures
  struct motif_group** motif_list = create_motif_list(N_loops,L_max); // indexed by n only, all motifs for that n
  struct motif_group*** motif_lib = create_motif_lib(N_loops,L_max); // indexed by n and L, only "lowest" rotation
  struct list* prev_motif_ids = create_list();

  // initialize with singleton motifs:
  for(int L=0;L<=L_max;L++) {
    struct motif* motif = create_motif(1);
    initialize_motif(motif,L); // one-loop motifs, setting l_0 = L;
    add_motif(motif_list[1],motif);
    append_motif(motif_lib[1][L],motif);
    list_append(prev_motif_ids,L);
  }
  
  // recursively build motif_lib for different lengths    
  for(int n=2;n<=N_loops;n++) { // create multi-junction loops 
    int new_size = 1 + L_max * pow(L_max+1,n-1); // 1 + largest possible id for n loops that add up to L_max or less
    printf("New motif size for n=%d: %d\n",n,new_size);
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
	  append_motif(motif_lib[n][L_tot],new_motif); //lexicographically lowest, only one added to lib
 	  list_append(motif_ids,new_id); // append the id of each motif
	  for(int j=0; j<n; j++) { // apart from repeats, motif with n loops has n rotations
	    struct motif* prev_motif = new_motif; // init. prev_motif to lexicographically lowest motif
	    new_motif = rotate_motif(new_motif,L_max); // first rotation
	    int new_id = get_motif_id(new_motif);
	    if(motif_empty(motif_list[n],new_id)) { // only add if it is not visited/stored yet
	      add_motif(motif_list[n],new_motif); // each rotation gets added to motif_list
	      list_append(motif_ids,new_id); // append the id of each rotation to motif_ids
	      set_next_motif(prev_motif,new_motif); // attach the rotation 
	    }
	  }
	}
      }
    }
    qsort(motif_ids->data,list_length(motif_ids),sizeof(int),compare);
    free_list(prev_motif_ids); // this should free memory allocated for motif_ids b/c of next line
    prev_motif_ids = motif_ids;
  }
  // final memory clean-up
  free_list(prev_motif_ids);
  free_motif_list(motif_list,N_loops);
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

void free_steptypes(struct steptypes* steps) {
  assert(steps);
  free(steps->steps);
  free(steps);
}

struct steptypes* motif_to_steptypes(struct motif* motif) {
  int num = get_motif_num(motif);
  int sum = get_motif_sum(motif);
  //printf("building motif n=%d,L=%d:",num,sum);
  //print_motif(motif,'\n');
  struct steptypes* steptypes = create_steptypes(2*num+sum);
  int count = 0;
  steptypes_set(steptypes,count,"BP");
  count++;
  for(int i=0; i<num; i++) {    
    for(int j=0; j<=motif->loop[i]; j++) {
      //printf("(i=%d,j=%d)",i,j);
      // loop of n nucleotides needs n+1 bb transforms, hence the "<=" in for loop
      steptypes_set(steptypes,count,"BB");      
      count++;
    }
    if(i < num-1) {
      //printf("[i=%d]",i);
      // add a base pair between each loop
      steptypes_set(steptypes,count,"BP");
      count++;
    }
  }
  //printf("\n");
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

struct steptypes* rotate_steptypes(struct steptypes* steptypes) {
  assert(steptypes);
  int length = steptypes_length(steptypes);
  struct steptypes* rotation = create_steptypes(length);
  // (a,b,c,...,y,z) -> (b,c,...,y,z,a)
  for(int i=1; i<length; i++) {
    steptypes_set(rotation,i-1,steptypes_get(steptypes,i));
  }
  steptypes_set(rotation,length-1,steptypes_get(steptypes,0));
  free_steptypes(steptypes); // free the incoming steptypes
  return rotation;
}

struct steptypes* steptypes_splice(struct steptypes* steptypes, int start, int end) {
  struct steptypes* substeps = create_steptypes(end-start+1); // positions are inclusive
  if(0 <= start && end < steptypes_length(steptypes)) {
    for(int i=start; i<=end; i++) {
      char* step = steptypes_get(steptypes,i);
      steptypes_set(substeps,i-start,step); // map positions to begin at 0
    }
  } else {
    printf("range error!");
    exit(0);
  }
  return substeps;
}

void print_steptypes(struct steptypes* steptypes) {
  for(int i=0; i<steptypes_length(steptypes); i++) {
    char* step = steptypes_get(steptypes,i);
    if(i == steptypes_length(steptypes)-1) {
      printf("%s",step);
    } else {
      printf("%s-",step);
    }
  }
  printf("\n");
}

void fprint_steptypes(FILE* F, struct steptypes* steptypes) {
  for(int i=0; i<steptypes_length(steptypes); i++) {
    char* step = steptypes_get(steptypes,i);
    if(i == steptypes_length(steptypes)-1) {
      fprintf(F,"%s",step);
    } else {
      fprintf(F,"%s-",step);
    }
  }
}
