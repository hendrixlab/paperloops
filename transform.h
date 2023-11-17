 /*
 * This library contains the definitions of structures and functions for transforms between 
 * adjacent nucleotides
 *
 * (C) David Hendrix david.hendrix@oregonstate.edu
 * In collaboration with Rhiju Das, Stanford
 * translated/derived from toyfold3 with permission: https://github.com/rhiju/toyfold3 (C) Rhiju Das
 */

#ifndef _TRANSFORM
#define _TRANSFORM

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "transform.h"
#include "motifs.h"
#include "rna_data.h"
#include "vectorfun.h"

/* STRUCTURES */

struct transform;

struct library;

/* FUNCTIONS */

// library tools 
struct library* create_library();
void free_library(struct library* transform_set);
int library_size(struct library* transform_set);
void resize_library(struct library* transform_set);
void remove_transform(struct library* transform_set, int idx);
struct transform* get_transform(struct library* transform_set, int idx);
struct transform* get_random_transform(struct library* lib);
void append_transform(struct library* transform_set, struct transform* T);

// frame tools 
struct frame* create_frame();
struct frame* get_frame(struct nuc* n);

// transform tools
struct transform* create_transform();
struct transform* compute_transform(struct frame* f1, struct frame* f2);
void set_euler_vector(struct transform* T);
double* get_euler_vector(struct transform* T);
double get_T6_value(struct transform* T, int i);
double* get_T6_std_dev(struct library* pts);
struct transform* apply_transform(struct transform* x1, struct transform* T);
struct transform* reverse_transform(struct transform* T);

// build libraries 
int valid_dinucleotide(struct nuc* n1, struct nuc* n2);
struct library* build_backbone_transform_library(struct rna_data* RNA);
struct library* build_basepair_transform_library(struct rna_data* RNA, struct bp_array* bp_array);

#endif
