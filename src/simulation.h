/*
 * This library contains the definitions of structures and functions for simulations and KDE
 * to perform basic statistics 
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 * In collaboration with Rhiju Das, Stanford
 * translated/derived from toyfold3 with permission: https://github.com/rhiju/toyfold3 (C) Rhiju Das
 */

#ifndef _SIMULATION
#define _SIMULATION

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "transform.h"
#include "motifs.h"
#include "vectorfun.h"
#include "rna_data.h"
#include "rna_io.h"

// this is just a static array of doubles for calculating averages/geometric means
struct array {
  double* data;
  int size;
};

// Data structures for collisions, self-avoiding polymers

// each monomer is a point in 3D space "r", and a numerical id number "num"
struct monomer {
  double* r;
  int num; // ID number avoids self-collisions
  struct monomer* next;
};

struct detector {
  struct monomer****  patch;
  int* bins;    // number of bins in each dimension
  int count;
  double* min;  // min positions in each dimension
  double* max;  // max position in each dimension
  double bin_size;
};

struct matrix {
  int* data;
  int rows;
  int cols;
  int transposed;
};

// subroutines for creating and using an array
struct array* create_array(int size);
void free_array(struct array* array);
int array_size(struct array* array);
void array_set(struct array* array, int idx, double val);
double array_get(struct array* array, int idx);

// basic statistics 
double get_sample_mean(double* samples, int N);
double get_sample_std_dev(double* samples, int N);
double get_geometric_mean(double val1, double val2);
struct array* get_samples_geometric_mean(struct array* samples1, struct array* samples2);

struct matrix* initialize_matrix(int rows, int cols);
void delete_matrix(struct matrix* M);
int get_element(struct matrix* M, int i, int j);
void set_element(struct matrix* M, int i, int j, int val);
void transpose_matrix(struct matrix* M);

struct monomer* create_monomer(double* x, int num);
void free_monomer(struct monomer* m);
struct detector* create_detector(struct trace* trace, double bin_size);
void free_detector(struct detector* detector);
int* get_patch_index(struct detector* detector, struct monomer* monomer);
double get_distance(struct monomer* a, struct monomer* b);
struct detector* build_detector(struct trace* trace, double bin_size);
struct detector* build_detector_with_helix(struct trace* trace, struct trace* helix, double bin_size);
void add_monomer_to_detector(struct detector* detector, struct monomer* new_mon);
void add_trace_to_detector(struct detector* detector, struct trace* trace);
int self_monomer_collision(struct detector* detector, struct monomer* mon_a, double d_min);
int monomer_collision(struct detector* detector, struct monomer* mon_a, double d_min);
int self_collision(struct detector* detector, struct trace* trace, int offset, double d_min);
int self_collision_helix(struct detector* detector, struct trace* trace, struct trace* helix, double d_min);
int trace_collision(struct detector* detector, struct trace* trace, double d_min);
double get_closest_approach(struct detector* detector, struct monomer* mon_a);

// simulation/sampling
struct trace* get_random_forward_trace(struct steptypes* steps, struct library* rna_lib);
struct trace* get_random_reverse_trace(struct steptypes* steps, struct library* rna_lib);
struct trace* build_random_forward_helix(struct steptypes* steps, struct trace* trace, struct library* rna_lib);
struct trace* build_random_reverse_helix(struct steptypes* steps, struct trace* trace, struct library* rna_lib);
struct trace* build_ideal_forward_helix(struct steptypes* steps, struct trace* trace, struct library* helix_lib);
struct trace* build_ideal_reverse_helix(struct steptypes* steps, struct trace* trace, struct library* helix_lib);
struct trace* build_forward_helix(struct steptypes* fwd_steps, struct trace* fwd_trace, struct library* rna_lib, struct library* helix_lib);
struct trace* build_reverse_helix(struct steptypes* rev_steps, struct trace* rev_trace, struct library* rna_lib, struct library* helix_lib);
struct trace* get_pts_forward(struct steptypes* steps, struct library* rna_lib, int size);
struct trace* get_pts_reverse(struct steptypes* steps, struct library* rna_lib, int size);
struct trace* get_pts_reverse2(struct steptypes* steps, struct library* rna_lib, int N_samples);

// KDE
struct array* get_kde_bandwidth(struct trace* pts);
struct array* get_kde_bandwidth_flat(struct array* pts, int N, double* h);
double gaussian_kernel(double x);
double get_kde_probability(struct frame* x, struct trace* pts, struct array* bandwidth);
double get_kde_probability_flat(int idx, struct array* pts1, struct array* pts2, double* bandwidth, int N2);
double get_kde_probability_sap(struct frame* x, struct trace* pts, struct array* bandwidth, struct matrix* M, int k);
double get_kde_probability_sap_flat(int idx, struct array* pts1, struct array* pts2, double* bandwidth, struct matrix* M, int N2);
double get_rna_closest_approach(struct rna_data* RNA, double bin_size);

// C_eff
double get_frame_dist_sq_flat(struct array* pts_fwd, int i, struct array* pts_rev, int j);
double get_frame_angle_sq_flat(struct array* pts_fwd, int i, struct array* pts_rev, int j);
double get_avg_dist_sq(struct array* pts_fwd, struct array* pts_rev, int N);
double get_avg_dist_sq_sap(struct array* pts_fwd, struct array* pts_rev, struct matrix* M, int N1, int N2);
double get_avg_angle_sq(struct array* pts_fwd, struct array* pts_rev, int N);
double get_avg_angle_sq_sap(struct array* pts_fwd, struct array* pts_rev, struct matrix* M, int N1, int N2);
struct trace* read_pts_file(char* pts_file);
void get_C_eff_samples(struct array* pts_fwd, struct array* pts_rev, double* h_rev, int N_samples, double* C_eff_samples);
int get_C_eff_samples_sap(struct array* pts_fwd, struct array* pts_rev, double* h_rev, struct matrix* M, int N_fwd, int N_rev, double* C_eff_samples);
struct array* get_C_eff_both_dir(struct steptypes* steps, struct library* rna_lib, int N_samples, double* C_eff_stats);
struct array* get_C_eff_both_dir_sap(struct motif* motif, struct steptypes* steps, struct library* rna_lib, int N_samples, double bin_size, double d_min, struct library* helix_lib, char* sap_report, double* C_eff_stats);
void fprint_loop_trace(FILE* OUT, int *idx, struct motif* motif, struct steptypes* fwd_steps, struct trace* fwd_trace, struct trace* fwd_helix, struct steptypes* rev_steps, struct trace* rev_trace, struct trace* rev_helix, double p_circ);
void print_loop_traces(FILE* OUT, int *idx, struct motif* motif, struct steptypes* steps, struct library* rna_lib, int N_samples, struct library* helix_lib, double p_circ_min);
void print_loop_traces_sap(FILE* OUT, int *idx, struct motif* motif, struct steptypes* steps, struct library* rna_lib, int N_samples, double bin_size, double d_min, struct library* helix_lib, char* sap_report, double p_circ_min);

#endif
