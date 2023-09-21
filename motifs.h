 /*
 * This library contains the definitions of structures and functions related to loop motifs
 * defining the lengths of each loop in a multiloop structure. e.g. [3,4,2,1]
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
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

struct motif;
struct motif_group;


/*
 * Tools for building/adding/removing to a motif
 */

struct motif* create_motif(int num);
void initialize_motif(struct motif* motif, int loop);
void set_motif_loop(struct motif* motif, int idx, int loop);
void set_motif_id(struct motif* motif, int id);
int get_motif_num(struct motif* motif);
int get_motif_loop(struct motif* motif, int idx);
int get_motif_id(struct motif* motif);
int get_motif_sum(struct motif* motif);
int print_motif(struct motif* motif, char c);
int get_new_id(int id, int new_num, int new_loop, int L_max);
int compute_motif_id(struct motif* motif, int L_max);
struct motif* expand_motif(struct motif* motif, int new_loop, int L_max);
struct motif* rotate_motif(struct motif* motif, int L_max);

/*
 * Tools for building/adding/removing elements from motif_group
 */

struct motif_group* create_motif_group(int size);
void free_motif_group(struct motif_group* motif_group);
int get_motif_group_size(struct motif_group* motif_group);
int get_motif_group_count(struct motif_group* motif_group);
struct motif* get_motif(struct motif_group* motif_group, int idx);
void add_motif(struct motif_group* motif_group, struct motif* m);
struct motif_group** build_motif_list(int N_loops, int L_max);

#endif


