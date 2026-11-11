 /*
 * This library contains the definitions of structures and functions related to loop motifs
 * defining the lengths of each loop in a multiloop structure. e.g. [3,4,2,1]
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 * In collaboration with Rhiju Das, Stanford
 * translated/derived from toyfold3 with permission: https://github.com/rhiju/toyfold3 (C) Rhiju Das
 */

#ifndef _MOTIFS
#define _MOTIFS

#include "motifs.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/*
 * data structures for the transform and for the library of transforms
 */

/*
 * data structures for the motifs and combinatorics
 */

// this is just a dynamic array of int values, used to store primes and divisors for various tasks
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

// a hash table structure to store motifs
struct hash_table {
  struct hash_node** array;
  int size;
  int total;
};

// this defines a single loop motif
struct motif {
  int* loop;
  int num;
  long int id; // a unique ID for each loop motif
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

/*
 * Tools for building/adding/removing to a motif
 */

int* create_rep_array(int len, int val);// list of integers using a dynamic array
struct list* create_list();
void free_list(struct list* list);
void list_append(struct list* list, long int n);
int list_length(struct list* list);
void print_list(struct list* list, char* title);

// hash table
int hash_function(struct hash_table* hash_table, long int key);
struct hash_table* create_hash_table(int array_size);
void free_hash_table(struct hash_table* hash_table);
void reset_hash_table(struct hash_table* hash_table);
void hash_table_set(struct hash_table* hash_table, long int key, void* value);
void* hash_table_get(struct hash_table* hash_table, long int key);
int hash_table_empty(struct hash_table* hash_table, long int key);

struct list* get_divisors(int n);
struct list*  get_primes(int n);
int mobius(int n, struct list*  primes);
int binomial(int n, int k);
int count_necklaces_fixed_sum(int n, int L);
int aperiodic_necklace_counts_fixed_length(int n, int L, struct list* divisors, struct list* primes);
int** create_necklace_counts(int N_loops, int L_max);
void free_necklace_counts(int** N, int N_loops);
int** necklace_counts_fixed_length(int N_loops, int L_max);

// motif tools 

struct motif* create_motif(int num);
void free_motif(struct motif* motif);
void initialize_motif(struct motif* motif, int loop);
void set_motif_loop(struct motif* motif, int idx, int loop);
void set_motif_id(struct motif* motif, long int id); 
int get_motif_num(struct motif* motif);
int get_motif_loop(struct motif* motif, int idx);
long int get_motif_id(struct motif* motif);
int get_motif_sum(struct motif* motif);
void print_motif(struct motif* motif, char c);
void fprint_motif(FILE* F, struct motif* motif, char c);
int sprint_motif(char* buffer, struct motif* motif, char c);
long int get_new_id(long int id, int new_num, int new_loop, int L_max);
long int compute_motif_id(struct motif* motif, int L_max);
struct motif* expand_motif(struct motif* motif, int new_loop, int L_max);
struct motif* rotate_motif(struct motif* motif, int L_max);
struct motif* get_next_motif(struct motif* motif);
struct motif* set_next_motif(struct motif* motif, struct motif* next);

/*
 * Tools for building/adding/removing elements from motif_group
 */

// motif_group - a list/collection of motifs
struct motif_group* create_motif_group(int size);
void free_motif_group(struct motif_group* motif_group);
int get_motif_group_size(struct motif_group* motif_group);
int get_motif_group_count(struct motif_group* motif_group);
struct motif* get_motif(struct motif_group* motif_group, int idx);
void set_motif(struct motif_group* motif_group, int idx, struct motif* motif);
void add_motif(struct motif_group* motif_group, struct motif* m);
void append_motif(struct motif_group* motif_group, struct motif* motif);

// motif_hash - a hash table for keeping track of motifs
struct hash_table** create_motif_hash(int N_loops, int L_max);
void free_motif_hash(struct hash_table** motif_hash, int N_loops);

// motif_lib - a library of motifs
struct motif_group*** create_motif_lib(int N_loops, int L_max);
void free_motif_lib(struct motif_group*** motif_lib, int N_loops, int L_max);
struct motif_group*** build_motif_lib(int N_min, int N_max, int L_max, int** N);

/***********************
 * Tools for Steptypes *
 ***********************/

struct steptypes* create_steptypes(int length);
void free_steptypes(struct steptypes* steps);
struct steptypes* clone_steptypes(struct steptypes* steps);
struct steptypes* motif_to_steptypes(struct motif* motif);
int motif_steptypes_length(struct motif* motif);
int steptypes_length(struct steptypes* steptypes);
int steptypes_bp_count(struct steptypes* steptypes);
char* steptypes_get(struct steptypes* steptypes, int idx);
void steptypes_set(struct steptypes* steptypes, int idx, char* step);
struct steptypes* rotate_steptypes(struct steptypes* steptypes);
struct steptypes* steptypes_splice(struct steptypes* steptypes, int start, int end);
struct steptypes* steptypes_concat(struct steptypes* steps1, struct steptypes* steps2);
void print_steptypes(struct steptypes* steptypes);
void fprint_steptypes(FILE* F, struct steptypes* steptypes);
int sprint_steptypes(char* buffer, struct steptypes* steptypes, char c);

#endif

