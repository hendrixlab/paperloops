 /*
 * This library contains the definitions of structures and functions for transforms between 
 * adjacent nucleotides
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 * In collaboration with Rhiju Das, Stanford
 * translated/derived from toyfold3 with permission: https://github.com/rhiju/toyfold3 (C) Rhiju Das
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "transform.h"
#include "vectorfun.h"

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
struct library* get_random_trace(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib);
struct library* get_pts_forward(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size);
struct library* get_pts_reverse(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size);

// KDE
struct array* get_kde_bandwidth(struct library* pts);
double get_kde_probability(struct transform* x, struct library* pts, double* bandwidth);

// C_eff
struct array* get_C_eff_samples(struct library* pts_fwd, struct library* pts_rev, struct array* h);
struct array* get_C_eff_both_dir(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size);
