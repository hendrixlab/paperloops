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
#include "rna_io.h"

// this is just an array of doubles
struct array {
  double* data;
  int size;
};

// this is creating a array, which is a dynamic array
struct array* create_array(int size) {
  struct array* array = malloc(sizeof(struct array));
  array->size = size;
  array->data = malloc(size * sizeof(double));
  return array;
}

void free_array(struct array* array) {
  free(array->data);
  free(array);
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
  append_transform(trace,x); // add the initial point

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

struct library* get_pts_forward(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int N_samples) {
  struct library* pts_fwd = create_library(); // create the library to store the points
  for(int i=0; i<N_samples; i++) {
    struct library* trace = get_random_trace(steps,bb_lib,bp_lib);
    struct transform* T_fwd = copy_transform(get_transform(trace,-1)); // get the last point in the trace
    double* EV = compute_euler_vector(get_rotation_matrix(T_fwd));
    set_euler_vector(T_fwd,EV);
    free_library(trace);
    append_transform(pts_fwd,T_fwd);
  }
  return pts_fwd;
}

struct library* get_pts_reverse(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int N_samples) {
  struct library* pts_fwd = get_pts_forward(steps,bb_lib,bp_lib,N_samples);
  struct library* pts_rev = create_library();
  for(int i=0; i<library_size(pts_fwd); i++) {
    struct transform* T_fwd = get_transform(pts_fwd,i); // T goes from start of steps to end for each replicate
    struct transform* T_rev = reverse_transform(T_fwd); // reverse it to go from end to start
    double* EV = compute_euler_vector(get_rotation_matrix(T_rev));
    set_euler_vector(T_rev,EV);
    append_transform(pts_rev,T_rev);
  }
  free_library(pts_fwd); // this should free the individual transforms in pts_fwd
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
    array_set(h, i, pow(4.0/(n*(d+2)),1.0/(d+4)) * std_dev[i] );
    //printf("bw calc: %d %f %f (d=%d,n=%d)\n",i,pow(4.0/(n*(d+2)),1.0/(d+4)), std_dev[i], d, n);
    //printf("%.6f\n",array_get(h,i));
  }  
  free(std_dev);
  return h;
}

double gaussian_kernel(double x) {
  return exp(-0.5*x*x)/sqrt(2.0*pi);
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

  double sum = 0.0;
  for(int i=0; i<n; i++) {
    double product = 1.0;
    struct transform* pts_T = get_transform(pts,i);
    for(int j=0; j<d; j++) {
      double x_val = get_T6_value(x,j);
      double pts_val = get_T6_value(pts_T,j);
      double h = array_get(bandwidth,j);
      double z = (x_val - pts_val)/h;
      product *= gaussian_kernel(z)/n;
      //printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n",i,j,x_val,pts_val,h,z,gaussian_kernel(z));
      //printf("%f\t%f\n",z,gaussian_kernel(z));
    }
    //printf("-> %.6e\n",product);
    //product /= n;
    sum += product;
  }
  //printf("sum -> %.6e\n",sum);
  sum /= det_bw;
  //printf("p = %.6e\n",sum);
  return sum;
}

/*
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
  double coeff = pow(2.0*pi,-d/2.0) / det_bw;
  
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
      printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\n",i,j,x_val,pts_val,h,val*val/(h*h),exponent);
    }
    sum += coeff*exp(-exponent/2.0);
  }
  return sum / (double)n;
}
*/

/*************************
 ** Collision Detection **
 *************************/

struct monomer {
  double* r;
  int num;
  struct monomer* next;
};

struct detector {
  struct monomer**** patch;
  int* bins;    // number of bins in each dimension
  int count;
  double* min;  // min positions in each dimension
  double* max;  // max position in each dimension
  double bin_size;
};

struct detector* create_detector(struct library* trace, double bin_size) {
  struct detector* detector = malloc(sizeof(struct detector));

  detector->bin_size = bin_size;
  detector->min = get_library_min_bounds(trace);
  detector->max = get_library_max_bounds(trace);
  detector->bins = malloc(3 * sizeof(int));
  // the number of bins in each dimension:
  for(int d=0; d<3; d++) {
    detector->bins[d] = (detector->max[d] - detector->min[d])/bin_size;
  }
  // a 3D array of "patches" defined as linked lists of monomers/atoms in that region
  detector->patch = malloc(detector->bins[0] * sizeof(struct monomer****));
  for(int i=0; i<detector->bins[0]; i++) {
    detector->patch[i] = malloc(detector->bins[1] * sizeof(struct monomer***));
    for(int j=0; j<detector->bins[1]; j++) {
      detector->patch[i][j] = malloc(detector->bins[2] * sizeof(struct monomer**));
      for(int k=0; k<detector->bins[2]; k++) {
	detector->patch[i][j][k] = NULL;
      }
    }
  }
  return detector;
}

int* get_patch_index(struct detector* detector, struct monomer* monomer) {
  int* index;
  for(int d=0; d<3; d++) {
    index[d] = (int) (monomer->r[d] - detector->min[d])/detector->bin_size;
  }
  return index;
}

double get_distance(struct monomer* a, struct monomer* b) {
  double* diff = subtract(a->r,b->r);
  double d = norm(diff);
  free_vector(diff);
  return d;
}

int collision(struct detector* detector, struct monomer* new_monomer) {
  double d_min = 0.5;
  int* idx = get_patch_index(detector,new_monomer);
  if(detector->patch[idx[0]][idx[1]][idx[2]]) {
    struct monomer* monomer = detector->patch[idx[0]][idx[1]][idx[2]];
    while(monomer) {
      double d = get_distance(monomer,new_monomer);
      if(d < d_min) {
	return 1;
      }
      monomer = monomer->next;
    }
  }
  free(idx);
  return 0;
}

void add_monomer(struct detector* detector, struct monomer* new_monomer) {
  int* idx = get_patch_index(detector,new_monomer);
  if(detector->patch[idx[0]][idx[1]][idx[2]]) {
    struct monomer* monomer = detector->patch[idx[0]][idx[1]][idx[2]];
    while(monomer->next) {
      monomer = monomer->next;
    }
    monomer->next = new_monomer;
  } else {
    detector->patch[idx[0]][idx[1]][idx[2]] = new_monomer;
  }
}

/*********************************
 * Effective Concentration C_eff * 
 *********************************/

struct array* get_C_eff_samples(struct library* pts_rev, struct library* pts_fwd, struct array* h_fwd) {
  int N = library_size(pts_rev);
  struct array* C_eff_samples = create_array(N);

  for(int i=0; i<N; i++) {
    struct transform* x = get_transform(pts_rev,i);

    /* testing block
    int d=6;
    printf("%.6f",get_T6_value(x,0));
    for(int j=1; j<d; j++) {
      printf("\t%.6f",get_T6_value(x,j));
    }
    printf("\n");
    // end test block */

    double p = get_kde_probability(x, pts_fwd, h_fwd);
    //printf("p = %.6e\n",p);
    double* EV = get_euler_vector(x);
    //printf("%.3f\t%.3f\t%.3f\n",EV[0],EV[1],EV[2]);
    double v = norm(get_euler_vector(x));
    //printf("%.6f\n",v);
    double w = 1.0;
    if(v != 0.0) {
      double sqrt_w = sin(v/2)/(v/2);
      w = sqrt_w*sqrt_w; // more efficient than pow()
      //printf("%.6f\n",w);
    }
    double C_eff = p/(w/(8.0*pi*pi))/(6.022e23/1.0e27);
    double den = (w/(8.0*pi*pi))/(6.022e23/1.0e27);
    //printf("C_eff = %.6f\n",C_eff);
    array_set(C_eff_samples,i,C_eff);
  }
  return C_eff_samples;
}


struct library* read_pts_file(char* pts_file) {
  int MAX_LINE_LENGTH = 1000;
  FILE* F = fopen(pts_file,"r");
  if(F == NULL) {
    printf("Error opening the file:\"%s\".\n Check permissions/file name.\n",pts_file);
    exit(1);
  }

  struct library* pts = create_library();
  char line[MAX_LINE_LENGTH];
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    struct split_text* row = split_line(line,"\t");
    if(get_split_text_cols(row) == 7) {
      if(strcmp(get_split_text_term(row,0),"idx") != 0) {
	struct transform* T = create_transform();
	double* t = create_vector();
	set_vector_elements(t,atof(get_split_text_term(row,1)),atof(get_split_text_term(row,2)),atof(get_split_text_term(row,3)));
	set_translation_vector(T,t);
	double* EV = create_vector();
	set_vector_elements(EV,atof(get_split_text_term(row,4)),atof(get_split_text_term(row,5)),atof(get_split_text_term(row,6)));
	//printf("test: %.3f %.3f %.3f\n",atof(get_split_text_term(row,4)),atof(get_split_text_term(row,5)),atof(get_split_text_term(row,6)));
	set_euler_vector(T,EV);
	append_transform(pts,T); 
      }
    }
  }
  return pts;
}

struct array* get_C_eff_both_dir(struct steptypes* steps, struct library* bb_lib, struct library* bp_lib, int N_samples) {
  int length = steptypes_length(steps);
  int halfway = length/2; // should be floored b/c ints

  //printf("length=%d, halfway=%d: [%d:%d] [%d:%d]\n",length,halfway,0,halfway-1,halfway,length-1);
  struct steptypes* steps_fwd = steptypes_splice(steps,0,halfway-1);
  struct library* pts_fwd =  get_pts_forward(steps_fwd,bb_lib,bp_lib,N_samples);
  
  struct steptypes* steps_rev = steptypes_splice(steps,halfway,length-1);
  struct library* pts_rev = get_pts_reverse(steps_rev,bb_lib,bp_lib,N_samples);

  /*
  //char* pts_fwd_file = "pts_fwd_file.txt";
  //char* pts_rev_file = "pts_rev_file.txt";
  char* pts_fwd_file = "../toyfold3/runs/sampled_pts_fwd.txt";
  char* pts_rev_file = "../toyfold3/runs/sampled_pts_rev.txt";

  struct library* pts_fwd = read_pts_file(pts_fwd_file);
  struct library* pts_rev = read_pts_file(pts_rev_file);
  
  /*
  int d = 6;  
  char pts_fwd_file[50] = "pts_fwd_file.txt";
  FILE* PTS_FWD_OUT = fopen(pts_fwd_file,"w");
  fprintf(PTS_FWD_OUT,"pl_x\tpl_y\tpl_z\tpl_vx\tpl_vy\tpl_vz\n");
  for(int i=0; i<library_size(pts_fwd); i++) {
    struct transform* T = get_transform(pts_fwd,i);
    set_euler_vector(T);
    double x = get_T6_value(T,0);
    fprintf(PTS_FWD_OUT,"%.6f",x);
    for(int j=1; j<d; j++) {
      double x = get_T6_value(T,j);
      fprintf(PTS_FWD_OUT,"\t%.6f",x);
    }
    fprintf(PTS_FWD_OUT,"\n");
  }
  fclose(PTS_FWD_OUT);
  
  char pts_rev_file[50] = "pts_rev_file.txt";
  FILE* PTS_REV_OUT = fopen(pts_rev_file,"w");
  fprintf(PTS_REV_OUT,"pl_x\tpl_y\tpl_z\tpl_vx\tpl_vy\tpl_vz\n");
  for(int i=0; i<library_size(pts_rev); i++) {
    struct transform* T = get_transform(pts_rev,i);
    set_euler_vector(T);
    double x = get_T6_value(T,0);
    fprintf(PTS_REV_OUT,"%.6f",x);
    for(int j=1; j<d; j++) {
      double x = get_T6_value(T,j);
      fprintf(PTS_REV_OUT,"\t%.6f",x);
    }
    fprintf(PTS_REV_OUT,"\n");
  }
  fclose(PTS_REV_OUT);
  */
  
  struct array* h_fwd = get_kde_bandwidth(pts_fwd);
  struct array* h_rev = get_kde_bandwidth(pts_rev);
  
  struct array* C_eff_fwd_samples = get_C_eff_samples(pts_rev,pts_fwd,h_fwd);
  struct array* C_eff_rev_samples = get_C_eff_samples(pts_fwd,pts_rev,h_rev);
  
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
  double C_eff_diff = 0.5 * fabs(log(C_eff_fwd_mean/C_eff_rev_mean));
  //printf("C_eff_f=%.3f C_eff_r=%.3f\n",C_eff_fwd_mean,C_eff_rev_mean);

  // geometric mean of the fwd and rev arithmetic means:
  double C_eff = get_geometric_mean(C_eff_fwd_mean,C_eff_rev_mean);
  // compute relative error as a Euclidian distance of statistical and systematic error "dimensions"
  double C_eff_relerr = sqrt(C_eff_fwd_relerr*C_eff_fwd_relerr + C_eff_rev_relerr*C_eff_rev_relerr + C_eff_diff*C_eff_diff);
  // compute error by scaling unitless relative error by the mean
  double C_eff_error = C_eff_relerr * C_eff;
  //printf("C_eff report: Geom-mean=%.3f Fwd-relerr=%.3f Rev-relerr=%.3f Syst-err=%.3f\n",C_eff,C_eff_fwd_relerr,C_eff_rev_relerr,C_eff_diff);
  //printf("C_eff_relerr: %.3f  ->> C_eff_err=%.3f\n",C_eff_relerr,C_eff_error);
  // set data to be returned 
  struct array* C_eff_stats = create_array(2);
  array_set(C_eff_stats,0,C_eff);
  array_set(C_eff_stats,1,C_eff_error);

  // memory managment
  free_steptypes(steps_fwd);
  free_steptypes(steps_rev);
  free_library(pts_fwd);
  free_library(pts_rev);
  free_array(h_fwd);
  free_array(h_rev);
  free_array(C_eff_fwd_samples);
  free_array(C_eff_rev_samples);
  
  return C_eff_stats;
}
