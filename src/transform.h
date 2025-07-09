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

struct frame;

struct trace;

struct transform;

struct library;

/* FUNCTIONS */

// frame tools
struct frame* create_frame();
void free_frame(struct frame* f);
struct frame* compute_frame(struct nuc* n);
struct frame* reverse_frame(struct frame* f);
struct frame* copy_frame(struct frame* f);
struct frame* average_frames(struct frame* f1, struct frame* f2);
double* get_frame_center(struct frame* f);
void set_frame_center(struct frame* f, double* ctr);
double** get_frame_orientation(struct frame* f);
void set_frame_orientation(struct frame* f, double** M);
void set_frame_euler_vector(struct frame* f, double* EV);
double* get_frame_euler_vector(struct frame* f);
int get_frame_num(struct frame* f);
void set_frame_num(struct frame* f, int num);
double get_T6_value(struct frame* f, int i);
void print_frame_T6(struct frame* f);
double* get_T6_std_dev(struct trace* pts);

// trace tools 
struct trace* create_trace(int length);
void free_trace(struct trace* trace);
int trace_length(struct trace* trace);
struct frame* get_frame(struct trace* trace, int idx);
void set_frame(struct trace* trace, int idx, struct frame* f);
void get_trace_min_bounds(struct trace* trace, double* min);
void get_trace_max_bounds(struct trace* trace, double* max);

// library tools 
struct library* create_library();
void free_library(struct library* lib);
int library_size(struct library* lib);
void resize_library(struct library* lib);
void remove_transform(struct library* lib, int idx);
struct transform* get_transform(struct library* lib, int idx);
struct transform* get_random_transform(struct library* lib);
void append_transform(struct library* lib, struct transform* T);
double* get_library_min_bounds(struct library* lib);
double* get_library_max_bounds(struct library* lib);

// transform tools
struct transform* create_transform(double val);
void free_transform(struct transform* T);
double* get_translation_vector(struct transform* T);
double get_transform_displacement(struct transform* T);
void set_translation_vector(struct transform* T, double* t);
double** get_rotation_matrix(struct transform* T);
void set_rotation_matrix(struct transform* T, double** R);
struct transform* copy_transform(struct transform* T);
struct transform* compute_transform(struct frame* f1, struct frame* f2);
//void set_euler_vector(struct transform* T, double* EV);
//double* get_euler_vector(struct transform* T);
//double get_T6_value(struct transform* T, int i);
//double* get_T6_std_dev(struct library* pts);
//void print_transform_T6(struct transform* T);
struct frame* apply_transform(struct transform* T, struct frame* x1);
struct transform* reverse_transform(struct transform* T);

// build libraries 
int valid_dinucleotide(struct nuc* n1, struct nuc* n2);
struct library* build_backbone_transform_library(struct rna_data* RNA);
struct library* build_basepair_transform_library(struct rna_data* RNA, struct bp_array* bp_array);

#endif
