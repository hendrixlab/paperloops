/*
 * This library contains the definitions of structures and functions for simulations and KDE
 * to perform basic statistics 
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 * In collaboration with Rhiju Das, Stanford
 * translated/derived from toyfold3 with permission: https://github.com/rhiju/toyfold3 (C) Rhiju Das
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "transform.h"
#include "motifs.h"
#include "vectorfun.h"
#include "rna_data.h"
#include "rna_io.h"

// this is just an array of doubles
struct array;

// useful subroutines for creating and using an array
struct array* create_array(int size);
void free_array(struct array* array);
int array_size(struct array* array);
void array_set(struct array* array, int idx, double val);
double array_get(struct array* array, int idx);

// basic statistics 
double get_sample_mean(struct array* samples);
double get_sample_std_dev(struct array* samples);
double get_geometric_mean(double val1, double val2);
struct array* get_samples_geometric_mean(struct array* samples1, struct array* samples2);

// simulation/sampling
struct trace* get_random_forward_trace(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib);
struct trace* get_random_reverse_trace(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib);
struct trace* build_random_forward_helices(struct steptypes* steps, struct trace* trace, struct library* bb_lib, struct library* bp_lib);
struct trace* build_random_reverse_helices(struct steptypes* steps, struct trace* trace, struct library* bb_lib, struct library* bp_lib);
struct trace* get_pts_forward(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size);
struct trace* get_pts_reverse(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size);
struct trace* get_pts_reverse2(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int N_samples);

// KDE
struct array* get_kde_bandwidth(struct library* pts);
double gaussian_kernel(double x);
double get_kde_probability(struct transform* x, struct library* pts, struct array* bandwidth);

// Collisions, self-avoiding polymers
struct monomer;
struct detector;

struct monomer* create_monomer(double* x, int num);
void free_monomer(struct monomer* m);
struct detector* create_detector(struct library* trace, double bin_size);
void free_detector(struct detector* detector);
int* get_patch_index(struct detector* detector, struct monomer* monomer);
double get_distance(struct monomer* a, struct monomer* b);
void add_monomer_to_detector(struct detector* detector, struct monomer* new_mon);
struct detector* add_trace_to_detector(struct detector* detector, struct trace* trace);
struct detector* build_detector(struct trace* trace, double bin_size);
int self_monomer_collision(struct detector* detector, struct monomer* mon_a, double d_min);
int monomer_collision(struct detector* detector, struct monomer* mon_a, double d_min);
int self_collision(struct detector* detector, struct library* trace, int offset, double d_min);
int trace_collision(struct detector* detector, struct trace* trace, double d_min);
double get_closest_approach(struct detector* detector, struct monomer* mon_a);
double get_rna_closest_approach(struct rna_data* RNA, double bin_size);

// C_eff
struct library* read_pts_file(char* pts_file);
struct array* get_C_eff_samples(struct library* pts_fwd, struct library* pts_rev, struct array* h);
struct array* get_C_eff_samples_sa(struct trace* pts_rev, struct trace* pts_fwd, struct array* h_fwd, int* rev_clear, int** rev_fwd_clear);
struct array* get_C_eff_both_dir(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size);
struct array* get_C_eff_both_dir_sa(struct motif* motif, struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int N_samples, double bin_size, double d_min, struct library* helix_bb_lib, struct library* helix_bp_lib, char* sa_report);

