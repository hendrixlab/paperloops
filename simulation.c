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
#include "vectorfun.h"

// this is just an array of doubles
struct array {
  double* data;
  int size;
};

// this is creating a array, which is a dynamic array
struct array* create_array(int size) {
  struct array* array = malloc(sizeof(array));
  array->size = size;
  array->data = malloc(size*sizeof(double));
  return array;
}

// return the current length of an array
int array_size(struct array* array) {
  return array->size;
}

// set a value to the array at index idx
void array_set(struct array* array, int idx, double val) {
  if(0 <= idx && idx < array_size(array)) {
    array->data[idx] = val;
  } else {
    printf("out of bounds in array_set()");
  }
}

// get a value from the array at index idx
double array_get(struct array* array, int idx) {
  if(0 <= idx && idx < array_size(array)) {
    return array->data[idx];
  } else {
    printf("out of bounds in array_get()");
  }
}

/********************
 * basic statistics *
 ********************/

double get_sample_mean(struct array* samples) {
  double sum = 0.0;
  int N = array_size(samples);
  for(int i=0; i<N; i++) {
    sum += array_get(samples,i);;
  }
  return sum / N;
}

double get_sample_std_dev(struct array* samples) {
  double sum = 0.0;
  double sum_squares = 0.0;
  int N = array_size(samples);
  
  for(int i=0; i<N; i++) {
     double sample = array_get(samples,i);
    sum += sample;
    sum_squares += sample*sample;
  }
  return sqrt( (sum_squares - sum*sum/N) / (N-1) );
}

double get_geometric_mean(double val1, double val2) {
  return sqrt(val1*val2);
}

struct array* get_samples_geometric_mean(struct array* samples1, struct array* samples2) {
  int N1 = array_size(samples1);
  int N2 = array_size(samples2);
  if(N1 == N2) {
    struct array* geometric_mean = create_array(N1);
    for(int i=0; i<N1; i++) {
      double val1 = array_get(samples1,i);;
      double val2 = array_get(samples2,i);;
      array_set(geometric_mean,i,get_geometric_mean(val1,val2));
    }
    return geometric_mean;
  } else {
    printf("Lengths do not match!\n");
    exit(0);
  }
}

/***********************
 * sampling/simulation *
 ***********************/

struct library* get_random_trace(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib) {
  struct library* trace = create_library(); // library to store the trace
  struct transform* x = create_transform(); // the initial point
  
  for(int i=0; i<steptypes_length(steps); i++) {
    char* step = steptypes_get(steps,i);
    struct transform* T;
    if(strcmp(step,"BB") == 0) {
      T = get_random_transform(bb_lib);
    } else if(strcmp(step,"BP") == 0) {
      T = get_random_transform(bp_lib);
    } else {
      printf("Unknown step found: %s\n",step);
      exit(1);
    }
    x = apply_transform(x,T);
    append_transform(trace,x);
  }
  return trace;
}

struct library* get_pts_forward(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size) {
  struct library* pts_fwd = create_library(); // create the library to store the points
  for(int i=0; i<size; i++) {
    struct library* trace = get_random_trace(steps,bb_lib,bp_lib);
    struct transform* x = get_transform(trace,-1); // get the last transform in the trace
    set_euler_vector(x);
    append_transform(pts_fwd,x);
  }
  return pts_fwd;
}

struct library* get_pts_reverse(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size) {
  struct library* pts_fwd = get_pts_forward(steps,bb_lib,bp_lib,size);
  struct library* pts_rev = create_library();
  for(int i=0; i<library_size(pts_fwd); i++) {
    struct transform* T = get_transform(pts_fwd,i);
    struct transform* T_rev = reverse_transform(T);
    set_euler_vector(T_rev);
    append_transform(pts_rev,T_rev);
  }
  return pts_rev;
}

/*************
 * KDE tools *
 *************/

struct array* get_kde_bandwidth(struct library* pts) {
  // calculate bandwidth using Silverman's rule
  int d = 6; // assuming 6D coordinates
  int n = library_size(pts); // number of points
  
  struct array* h = create_array(d);
  double* std_dev = get_T6_std_dev(pts);
  
  for(int i=0; i<d; i++) {
    // the h[i] terms correspond to diagonal terms sqrt(H[i,i]) of a bandwidth matrix
    array_set(h, i, pow(4/(n*(d+2)),1/(n+4)) * std_dev[i] );
  }  
  free(std_dev);
  return h;
}

double get_kde_probability(struct transform* x, struct library* pts, struct array* bandwidth) {
  // evaluate the probability f(x) of point x using the data in pts
  // essentially follows https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation
  // using standard multivariate normal kernels
  int d = 6;
  int n = library_size(pts);
  
  // bandwidth is an array, terms of diagonal H matrix
  // determinant is the product of the diagonal terms
  double det_bw = 1.0;
  for(int j=0; j<d; j++) {
    det_bw *= array_get(bandwidth,j);
  }
  double coeff = pow(2.0*pi,-d/2) / det_bw;
  
  double sum = 0.0;
  for(int i=0; i<n; i++) {
    struct transform* pts_T = get_transform(pts,i);
    //printf("pts_T_%d=",i);
    //print_transform_T6(pts_T);
    double exponent = 0.0; // will store the sum resulting from (x-x')H(x-x')
    for(int j=0; j<d; j++) {
      double x_val = get_T6_value(x,j);
      double pts_val = get_T6_value(pts_T,j);
      double val = x_val - pts_val;
      double h = array_get(bandwidth,j);
      exponent += val*val/(h*h);
    }
    sum += coeff*exp(-exponent/2.0);
  }
  return sum / (double)n;
}

/*********************************
 * Effective Concentration C_eff * 
 *********************************/

struct array* compute_C_eff_samples(struct library* pts_fwd, struct library* pts_rev, struct array* h) {
  int N = library_size(pts_rev);
  struct array* C_eff_samples = create_array(N);
  
  for(int i=0; i<N; i++) {
    struct transform* x = get_transform(pts_rev,i);
    double p = get_kde_probability(x, pts_fwd, h);
    //printf("p=%.6e\n",p);
    double v = norm(get_euler_vector(x));
    double w = 1.0;
    if(v != 0.0) {
      double sqrt_w = sin(v/2)/(v/2);
      w = sqrt_w*sqrt_w; // more efficient than pow()
    }
    double C_eff = p/(w/(8*pi*pi))/(6.022e23/1e27);
    array_set(C_eff_samples,i,C_eff);
  }
  return C_eff_samples;
}

struct array* get_C_eff_both_dir(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int size) {
  int length = steptypes_length(steps);
  int halfway = length/2; // should be floored b/c ints

  //printf("length=%d, halfway=%d: [%d:%d] [%d:%d]\n",length,halfway,0,halfway-1,halfway,length-1);
  struct steptypes* steps_fwd = steptypes_splice(steps,0,halfway-1);
  struct library* pts_fwd =  get_pts_forward(steps_fwd,bb_lib,bp_lib,size);
  
  struct steptypes* steps_rev = steptypes_splice(steps,halfway,length-1);
  struct library* pts_rev = get_pts_reverse(steps_rev,bb_lib,bp_lib,size);
  
  struct array* h_fwd = get_kde_bandwidth(pts_fwd);
  struct array* h_rev = get_kde_bandwidth(pts_rev);

  struct array* C_eff_fwd_samples = compute_C_eff_samples(pts_fwd,pts_rev,h_fwd);
  struct array* C_eff_rev_samples = compute_C_eff_samples(pts_rev,pts_fwd,h_rev);
  
  double C_eff_fwd_mean = get_sample_mean(C_eff_fwd_samples);
  double C_eff_fwd_std_dev = get_sample_std_dev(C_eff_fwd_samples);
  double C_eff_fwd_error = C_eff_fwd_std_dev / sqrt(array_size(C_eff_fwd_samples));
  double C_eff_rev_mean = get_sample_mean(C_eff_rev_samples);
  double C_eff_rev_std_dev = get_sample_std_dev(C_eff_rev_samples);
  double C_eff_rev_error = C_eff_rev_std_dev / sqrt(array_size(C_eff_rev_samples));
  
  // statistical error
  double C_eff_fwd_relerr = C_eff_fwd_error / C_eff_fwd_mean;
  double C_eff_rev_relerr = C_eff_rev_error / C_eff_rev_mean;
  // systematic error
  double C_eff_diff = 0.5 * abs(log(C_eff_fwd_mean/C_eff_rev_mean));

  // geometric mean of the fwd and rev arithmetic means:
  double C_eff = get_geometric_mean(C_eff_fwd_mean,C_eff_rev_mean);
  // compute relative error as a Euclidian distance of statistical and systematic error "dimensions"
  double C_eff_relerr = sqrt(C_eff_fwd_relerr*C_eff_fwd_relerr + C_eff_rev_relerr*C_eff_rev_relerr + C_eff_diff*C_eff_diff);
  // compute error by scaling unitless relative error by the mean
  double C_eff_error = C_eff_relerr * C_eff;

  struct array* C_eff_stats = create_array(2);
  array_set(C_eff_stats,0,C_eff);
  array_set(C_eff_stats,1,C_eff_error);
  
  return C_eff_stats;
}
