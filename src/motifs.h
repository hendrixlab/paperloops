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

struct list;
struct hash_node;
struct hash_table;
struct motif;
struct motif_group;
struct motif_library;
struct steptypes;

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
struct hash_table* hash_table_create(int array_size);
void hash_table_free(struct hash_table* hash_table);
void hash_table_reset(struct hash_table* hash_table);
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
long int get_new_id(long int id, int new_num, int new_loop, int L_max);
long int compute_motif_id(struct motif* motif, int L_max);
struct motif* expand_motif(struct motif* motif, int new_loop, int L_max);
struct motif* rotate_motif(struct motif* motif, int L_max);
struct motif* get_next_motif(struct motif* motif);

/*
 * Tools for building/adding/removing elements from motif_group
 */

struct motif_group* create_motif_group(int size);
void free_motif_group(struct motif_group* motif_group);
int get_motif_group_size(struct motif_group* motif_group);
int get_motif_group_count(struct motif_group* motif_group);
struct motif* get_motif(struct motif_group* motif_group, int idx);
void add_motif(struct motif_group* motif_group, struct motif* m);
void append_motif(struct motif_group* motif_group, struct motif* motif);
struct motif_group*** create_motif_lib(int N_loops, int L_max);
void free_motif_lib(struct motif_group*** motif_lib, int N_loops, int L_max);
struct hash_table** create_motif_list(int N_loops, int L_max);
//void free_motif_list(struct motif_group** motif_list, int N_loops);
void free_motif_list(struct hash_table** motif_list, int N_loops);
struct motif_group*** build_motif_lib(int N_loops, int L_max, int** N);

/***********************
 * Tools for Steptypes *
 ***********************/

struct steptypes* create_steptypes(int length);
void free_steptypes(struct steptypes* steps);
struct steptypes* motif_to_steptypes(struct motif* motif);
int steptypes_length(struct steptypes* steptypes);
char* steptypes_get(struct steptypes* steptypes, int idx);
void steptypes_set(struct steptypes* steptypes, int idx, char* step);
struct steptypes* rotate_steptypes(struct steptypes* steptypes);
struct steptypes* steptypes_splice(struct steptypes* steptypes, int start, int end);
void print_steptypes(struct steptypes* steptypes);
void fprint_steptypes(FILE* F, struct steptypes* steptypes);

#endif

