 /*
 * This library contains the definitions of structures and functions for transforms between 
 * adjacent nucleotides
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
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

struct library* create_library();
void free_library(struct library* transform_set);
int library_size(struct library* transform_set);
void resize_library(struct library* transform_set);
void remove_transform(struct library* transform_set, int idx);
struct transform* get_transform(struct library* transform_set, int idx);
void add_transform(struct library* transform_set, struct transform* T);
struct frame* create_frame();
struct frame* get_frame(struct nuc* n);
struct transform* compute_transform(struct frame* f1, struct frame* f2);
int valid_dinucleotide(struct nuc* n1, struct nuc* n2);
struct library* build_backbone_transform_library(struct rna_data* RNA);
struct library* build_basepair_transform_library(struct rna_data* RNA, struct bp_array* bp_array);
struct transform* apply_transform(struct transform* x1, struct transform* T);
struct library* get_random_trace(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib);
struct library* get_pts_forward(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size);


#endif
