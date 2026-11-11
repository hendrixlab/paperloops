
/*
 * This library contains the definitions of structures and functions for simulations and KDE
 * to perform basic statistics 
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 * In collaboration with Rhiju Das, Stanford
 * translated/derived from toyfold3 with permission: https://github.com/rhiju/toyfold3 (C) Rhiju Das
 */
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#include "simulation.h"
#include "transform.h"
#include "motifs.h"
#include "vectorfun.h"
#include "rna_data.h"
#include "rna_io.h"

// this is creating static array of doubles 
struct array* create_array(int size) {
  struct array* array = malloc(sizeof(struct array));
  array->size = size;
  array->data = malloc(size * sizeof(double));
  return array;
}

// free the array of doubles 
void free_array(struct array* array) {
  assert(array);
  free(array->data);
  free(array);
}

// return the current length of an array
int array_size(struct array* array) {
  return array->size;
}

// set a value to the array at index idx
void array_set(struct array* array, int idx, double val) {
  assert(array);
  if(0 <= idx && idx < array_size(array)) {
    array->data[idx] = val;
  } else {
    printf("out of bounds in array_set()");
  }
}

// get a value from the array at index idx
double array_get(struct array* array, int idx) {
  assert(array);
  if(0 <= idx && idx < array_size(array)) {
    return array->data[idx];
  } else {
    printf("out of bounds in array_get()");
  }
}

/********************
 * basic statistics *
 ********************/

double get_mean(double val1, double val2) {
  return (val1+val2)/2.0;
}

double get_sample_mean(double* samples, int N) {
  assert(samples);
  double sum = 0.0;
  for(int i=0; i<N; i++) {
    sum += samples[i];
  }
  return sum / N;
}

double get_sample_std_dev(double* samples, int N) {
  assert(samples);
  double sum = 0.0;
  double sum_squares = 0.0;
  
  for(int i=0; i<N; i++) {
    double sample = samples[i];
    sum += sample;
    sum_squares += sample*sample;
  }
  return sqrt( (sum_squares - sum*sum/N) / (N-1) );
}

double get_geometric_mean(double val1, double val2) {
  return sqrt(val1*val2);
}

struct array* get_samples_geometric_mean(struct array* samples1, struct array* samples2) {
  assert(samples1);
  assert(samples2);
  int N1 = array_size(samples1);
  int N2 = array_size(samples2);
  if(N1 == N2) {
    struct array* geometric_mean = create_array(N1);
    for(int i=0; i<N1; i++) {
      double val1 = array_get(samples1,i);
      double val2 = array_get(samples2,i);
      array_set(geometric_mean,i,get_geometric_mean(val1,val2));
    }
    return geometric_mean;
  } else {
    printf("Error: Lengths do not match in get_samples_geometric_mean\n");
    exit(0);
  }
}

/*************************
 ** Collision Detection **
 *************************/

struct matrix* initialize_matrix(int rows, int cols) {
  struct matrix* M = malloc(sizeof(struct matrix));
  M->data = malloc(rows * cols * sizeof(int));
  M->rows = rows;
  M->cols = cols;
  M->transposed = 0;
  for(int i=0; i<rows*cols; i++) {
    M->data[i] = 1;
  }
  return M;
}

void delete_matrix(struct matrix* M) {
  free(M->data);
  free(M);
}

int get_element(struct matrix* M, int i, int j) {
  if(M->transposed) {
    return M->data[(j * M->cols) + i];
  }
  return M->data[(i * M->cols) + j];
}

void set_element(struct matrix* M, int i, int j, int val) {
  if(M->transposed) {
    M->data[(j * M->cols) + i] = val;
  } else {
    M->data[(i * M->cols) + j] = val;
  }
}

void transpose_matrix(struct matrix* M) {
  if(M->transposed) {
    M->transposed = 0;
  }
  M->transposed = 1;
}

/***********************************
 * Monomers, collisions, detectors *
 ***********************************/

struct monomer* create_monomer(double* x, int num) {
  struct monomer* monomer = malloc(sizeof(struct monomer));
  monomer->r = x;
  monomer->num = num;
  monomer->next = NULL;
  return monomer;
}

void free_monomer(struct monomer* m) {
  assert(m);  
  free(m);
}

struct detector* create_detector(struct trace* trace, double bin_size) {
  assert(trace);  
  struct detector* detector = malloc(sizeof(struct detector));

  detector->count = 0;
  detector->bin_size = bin_size;
  double* min = get_trace_min_bounds(trace,NULL);
  double* max = get_trace_max_bounds(trace,NULL);
  detector->min = min;
  detector->max = max;

  detector->bins = malloc(3 * sizeof(int));
  // set the number of bins in each dimension:
  for(int d=0; d<3; d++) {
    detector->bins[d] = ceil((detector->max[d] - detector->min[d])/bin_size);
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

struct detector* create_detector_with_helix(struct trace* trace, struct trace* helix, double bin_size) {
  assert(trace);
  assert(trace);
  struct detector* detector = malloc(sizeof(struct detector));

  detector->count = 0;
  detector->bin_size = bin_size;
  double* min = get_trace_min_bounds(trace,NULL); // define the bounds to include both trace and helix
  double* max = get_trace_max_bounds(trace,NULL);
  min = get_trace_min_bounds(helix,min);
  max = get_trace_max_bounds(helix,max);
  detector->min = min;
  detector->max = max;

  detector->bins = malloc(3 * sizeof(int));
  // set the number of bins in each dimension:
  for(int d=0; d<3; d++) {  
    detector->bins[d] = ceil((detector->max[d] - detector->min[d])/bin_size);
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

void free_detector(struct detector* detector) {  
  assert(detector);
  for(int i=0; i<detector->bins[0]; i++) {
    for(int j=0; j<detector->bins[1]; j++) {
      for(int k=0; k<detector->bins[2]; k++) {
	if(detector->patch[i][j][k]) {
	  struct monomer* m = detector->patch[i][j][k];
	  while(m->next) {
	    struct monomer* n = m->next;
	    free_monomer(m);
	    m = n;
	  }
	  free_monomer(m);
	}
      }
      free(detector->patch[i][j]);
    }
    free(detector->patch[i]);
  }
  free_vector(detector->min);
  free_vector(detector->max);
  free(detector->bins);
  free(detector->patch);
  free(detector);
}

// detector tools

int* get_patch_index(struct detector* detector, struct monomer* monomer) {
  int* index = malloc(3*sizeof(int));
  for(int d=0; d<3; d++) {
    index[d] = (int)((monomer->r[d] - detector->min[d])/detector->bin_size);
  }
  return index;
}

double get_distance(struct monomer* a, struct monomer* b) {
  double* diff = subtract(a->r,b->r);
  double d = norm(diff);
  free_vector(diff);
  return d;
}

int check_index(struct detector* detector, int* idx) {
  assert(detector);
  for(int d=0; d<3; d++) {
    if((idx[d] < 0)||(detector->bins[d] <= idx[d])) {
      printf("index out of bound for dimension %d: %d vs %d",d,idx[d],detector->bins[d]);
      exit(1);
    }
  }
  if(detector->patch[idx[0]][idx[1]][idx[2]]) {
    return 1;
  } else {
    return 0;
  }
}

// detector building

struct detector* build_detector(struct trace* trace, double bin_size) {
  assert(trace);
  struct detector* detector = create_detector(trace,bin_size);
  add_trace_to_detector(detector,trace);
  return detector;
}

struct detector* build_detector_with_helix(struct trace* trace, struct trace* helix, double bin_size) {
  assert(trace);
  assert(helix);
  struct detector* detector = create_detector_with_helix(trace,helix,bin_size);
  add_trace_to_detector(detector,trace);
  add_trace_to_detector(detector,helix);  
  return detector;
}

void add_monomer_to_detector(struct detector* detector, struct monomer* new_mon) {
  assert(detector);
  assert(new_mon);

  int* idx = get_patch_index(detector,new_mon);
  
  if(check_index(detector,idx)) {
    struct monomer* monomer = detector->patch[idx[0]][idx[1]][idx[2]];
    while(monomer->next) {
      monomer = monomer->next;
    }
    monomer->next = new_mon;
  } else {
    detector->patch[idx[0]][idx[1]][idx[2]] = new_mon;
  }
  detector->count++;
  free(idx);
}

void add_trace_to_detector(struct detector* detector, struct trace* trace) {
  assert(detector);
  assert(trace);
  int count = detector->count;
  for(int i=0; i<trace_length(trace); i++) {
    struct frame* f = get_frame(trace,i);
    double* x = get_frame_center(f);
    // update monomer index to distinguish traces from each other
    struct monomer* m = create_monomer(x,count+i);
    add_monomer_to_detector(detector,m);
  }
}

// actual collisions:

int monomer_collision(struct detector* detector, struct monomer* mon_a, double d_min) {
  // checks if a monomer collides with anything stored in the detector
  // no assumption on 
  assert(detector);
  assert(mon_a);
  //double d_min = 3;
  int* ctr_idx = get_patch_index(detector,mon_a);
  int* idx = malloc(3*sizeof(int));
  for(int shift0=-1; shift0<=1; shift0++) {
    for(int shift1=-1; shift1<=1; shift1++) {
      for(int shift2=-1; shift2<=1; shift2++) {
	idx[0] = ctr_idx[0] + shift0;
	idx[1] = ctr_idx[1] + shift1;
	idx[2] = ctr_idx[2] + shift2;
	if((0 <= idx[0])&&(idx[0] < detector->bins[0])&&
	   (0 <= idx[1])&&(idx[1] < detector->bins[1])&&
	   (0 <= idx[2])&&(idx[2] < detector->bins[2])) {
	  struct monomer* mon_b = detector->patch[idx[0]][idx[1]][idx[2]];
	  while(mon_b) {
	    // if both num=0 it is the origin
	    if(mon_a->num != 0 || mon_b->num != 0) {
	      double d = get_distance(mon_a,mon_b);
	      if(d < d_min) {
		free(idx);
		free(ctr_idx);
		return 1;
	      }
	    }
	    mon_b = mon_b->next;
	  }
	}
      }
    }
  }
  free(idx);
  free(ctr_idx);
  return 0;
}

int self_monomer_collision(struct detector* detector, struct monomer* mon_a, double d_min) {
  assert(detector);
  assert(mon_a);
  //double d_min = 3.0;
  int* ctr_idx = get_patch_index(detector,mon_a);
  int* idx = malloc(3*sizeof(int));
  for(int shift0=-1; shift0<=1; shift0++) {
    for(int shift1=-1; shift1<=1; shift1++) {
      for(int shift2=-1; shift2<=1; shift2++) {
	idx[0] = ctr_idx[0] + shift0;
	idx[1] = ctr_idx[1] + shift1;
	idx[2] = ctr_idx[2] + shift2;
	if((0 <= idx[0])&&(idx[0] < detector->bins[0])&&
	   (0 <= idx[1])&&(idx[1] < detector->bins[1])&&
	   (0 <= idx[2])&&(idx[2] < detector->bins[2])) {	  
	  struct monomer* mon_b = detector->patch[idx[0]][idx[1]][idx[2]];
	  while(mon_b) {
	    if(mon_a->num != mon_b->num) { // this if-statement for self collisions, avoid comparing to itself
	      double d = get_distance(mon_a,mon_b);
	      if(d < d_min) {
		free(idx);
		free(ctr_idx);
		return 1;
	      }
	    }
	    mon_b = mon_b->next;
	  }
	}
      }
    }
  }
  free(idx);
  free(ctr_idx);
  return 0;
}

int self_collision(struct detector* detector, struct trace* trace, int offset, double d_min) {
  // offset shifts the ID of the incoming monomer if helix and loops use same numbering system
  assert(detector);
  assert(trace);  
  // It is assumed that detector was built with the same trace 
  // functionally a boolean -- does this trace have a collision with itself?
  // int offset is the shift to the monomer IDs, used for helix
  for(int i=0; i<trace_length(trace); i++) {
    struct frame* f = get_frame(trace,i);
    double* x = get_frame_center(f);
    struct monomer* m = create_monomer(x,offset+i);
    if(self_monomer_collision(detector,m,d_min)) {
      // collides with current polymer
      free_monomer(m);
      return 1;
    } else {
      free_monomer(m);
    }
  }
  return 0;
}

int self_collision_helix(struct detector* detector, struct trace* trace, struct trace* helix, double d_min) {
  if(self_collision(detector,trace,0,d_min)) {
    return 1;
  } else if(self_collision(detector,helix,trace_length(trace),d_min)) {
    return 1;
  }
  return 0;  
}

int trace_collision(struct detector* detector, struct trace* trace, double d_min) {
  // for the pre-built detector test if the given trace collides
  assert(detector);
  assert(trace);
  for(int i=0; i<trace_length(trace); i++) {
    struct frame* f = get_frame(trace,i);
    double* x = get_frame_center(f);
    struct monomer* m = create_monomer(x,i);
    if(monomer_collision(detector,m,d_min)) {
      free_monomer(m);
      return 1;
    } else {
      free_monomer(m);
    }
  }
  return 0;  
}

int trace_collision_helix(struct detector* detector, struct trace* trace, struct trace* helix, double d_min) {
  if(trace_collision(detector,trace,d_min)) {
    return 1;
  } else if(trace_collision(detector,helix,d_min)) {
    return 1;
  }
  return 0;  
}

double get_closest_approach(struct detector* detector, struct monomer* mon_a) {
  assert(detector);
  assert(mon_a);
  
  double d_min = 1000000.0;
  int* ctr_idx = get_patch_index(detector,mon_a);
  int* idx = malloc(3*sizeof(int));
  for(int shift0=-1; shift0<=1; shift0++) {
    for(int shift1=-1; shift1<=1; shift1++) {
      for(int shift2=-1; shift2<=1; shift2++) {
	idx[0] = ctr_idx[0] + shift0;
	idx[1] = ctr_idx[1] + shift1;
	idx[2] = ctr_idx[2] + shift2;
	if((0 <= idx[0])&&(idx[0] < detector->bins[0])&&
	   (0 <= idx[1])&&(idx[1] < detector->bins[1])&&
	   (0 <= idx[2])&&(idx[2] < detector->bins[2])) {
	  struct monomer* mon_b = detector->patch[idx[0]][idx[1]][idx[2]];
	  while(mon_b) {
	    if(mon_a->num != mon_b->num) {
	      double d = get_distance(mon_a,mon_b);
	      if(d < d_min) {
		d_min = d;
		//printf("updating d_min to %.3f %d vs %d\n",d_min,mon_a->num,mon_b->num);
	      }
	    }
	    mon_b = mon_b->next;
	  }
	}
      }
    }
  }
  free(idx);
  free(ctr_idx);
  return d_min;
}

/***********************
 * sampling/simulation *
 ***********************/

struct trace* get_random_forward_trace(struct steptypes* steps, struct library* rna_lib) {
  assert(steps);
  assert(rna_lib);
  assert(rna_lib->bb);
  assert(rna_lib->bp);
  int length = steptypes_length(steps);
  struct trace* trace = create_trace(length+1); // structure to store the trace +1 for initial point stored at 0
  struct frame* x = create_frame(1.0); // the initial point
  set_frame(trace,0,x); // add the initial point

  for(int i=1; i<=steptypes_length(steps); i++) {
    char* step = steptypes_get(steps,i-1);
    struct transform* T;
    if(strcmp(step,"BB") == 0) {
      T = get_random_transform(rna_lib->bb);
    } else if(strcmp(step,"BP") == 0) {
      T = get_random_transform(rna_lib->bp);
    } else {
      printf("Unknown step found: %s\n",step);
      exit(1);
    }
    x = apply_transform(T,x);
    set_frame(trace,i,x);
  }
  return trace;
}

struct trace* get_random_reverse_trace(struct steptypes* steps, struct library* rna_lib) {
  assert(steps);
  assert(rna_lib);
  assert(rna_lib->bb);
  assert(rna_lib->bp);
  // input "steps" are reversed within this function, presumed input in forward direction
  int length = steptypes_length(steps);
  struct trace* trace = create_trace(length+1); // structure to store the trace
  struct frame* x = create_frame(1.0); // the initial point
  set_frame(trace,0,x); // add the initial point

  for(int i=steptypes_length(steps)-1; i>=0; i--) {
    char* step = steptypes_get(steps,i);
    struct transform* T;
    if(strcmp(step,"BB") == 0) {
      T = get_random_transform(rna_lib->bb);
    } else if(strcmp(step,"BP") == 0) {
      T = get_random_transform(rna_lib->bp);
    } else {
      printf("Unknown step found: %s\n",step);
      exit(1);
    }
    // reverse each individual transform in reverse order
    struct transform* Tr = reverse_transform(T);
    x = apply_transform(Tr,x);
    free_transform(Tr);
    set_frame(trace,steptypes_length(steps)-i,x); // these are stored in reverse order so that connecting point is last (0,0,0) first
  }
  return trace;
}

struct trace* build_random_forward_helix(struct steptypes* steps, struct trace* trace, struct library* rna_lib) {
  /*
                 Tbp2     
              f3  ->  f4
              ^       |
          Tf2 |  Tbp1 V
              f1  ->  f2
          Tf1 ^       | Tr
              |       V
      fi-1 -> fi  ->  fi+1 -> ...
        T[i-1]   T[i]    T[i+1]
  */
  assert(steps);
  assert(trace);
  assert(rna_lib);
  assert(rna_lib->bb);
  assert(rna_lib->bp);
  int Nbp = steptypes_bp_count(steps);
  //struct trace* helix = create_trace(Nbp*4);
  struct trace* helix = create_trace(Nbp*6);
  Nbp = 0;
  for(int i=0; i<steptypes_length(steps); i++) {
    char* step = steptypes_get(steps,i);
    if(strcmp(step,"BP") == 0) {
      struct frame* f = get_frame(trace,i);
      struct transform* Tf1 = get_random_transform(rna_lib->bb);
      struct frame* f1 = apply_transform(Tf1,f);
      struct transform* Tbp1 = get_random_transform(rna_lib->bp);
      struct frame* f2 = apply_transform(Tbp1,f1);
      struct transform* Tf2 = get_random_transform(rna_lib->bb);
      struct frame* f3 = apply_transform(Tf2,f1);
      struct transform* Tbp2 = get_random_transform(rna_lib->bp);
      struct frame* f4 = apply_transform(Tbp2,f3);
      struct frame* f12 = average_frames(f1,f2);
      struct frame* f34 = average_frames(f3,f4);
      set_frame(helix,6*Nbp,f1);
      set_frame(helix,6*Nbp+1,f2);
      set_frame(helix,6*Nbp+2,f12);
      set_frame(helix,6*Nbp+3,f3);
      set_frame(helix,6*Nbp+4,f4);
      set_frame(helix,6*Nbp+5,f34);
      Nbp++;
    }
  }
  return helix;
}

struct trace* build_random_reverse_helix(struct steptypes* steps, struct trace* trace, struct library* rna_lib) {
  /*
                 Tbp2     
              f4  <-  f3
              |       ^
          Tf2 V  Tbp1 |
              f2  <-  f1
          Tf1 |       ^ Tr
              V       |
      fi-1 <- fi  <-  fi+1 <- ...
        T[i-1]   T[i]    T[i+1]
  */
  assert(steps);
  assert(trace);
  assert(rna_lib);
  assert(rna_lib->bb);
  assert(rna_lib->bp);
  int Nbp = steptypes_bp_count(steps);
  struct trace* helix = create_trace(Nbp*6);
  Nbp = 0;
  if(steptypes_length(steps)+1 != trace_length(trace)) {
    printf("error in lengths in build_random_reverse_helix: %d vs %d\n",steptypes_length(steps),trace_length(trace));
    exit(1);
  }
  for(int i=steptypes_length(steps)-1; i>=0; i--) {
    char* step = steptypes_get(steps,i);
    if(strcmp(step,"BP") == 0) {
      struct frame* f = get_frame(trace,i-(steptypes_length(steps)-1)); // there should be 1 more in trace vs steps
      struct transform* Tf1 = get_random_transform(rna_lib->bb);
      struct transform* Tr1 = reverse_transform(Tf1);
      struct frame* f1 = apply_transform(Tr1,f);
      free_transform(Tr1);
      struct transform* Tbp1 = get_random_transform(rna_lib->bp);
      struct transform* Trbp1 = reverse_transform(Tbp1);
      free_transform(Trbp1);
      struct frame* f2 = apply_transform(Trbp1,f1);
      struct transform* Tf2 = get_random_transform(rna_lib->bb);
      struct transform* Tr2 = reverse_transform(Tf2);
      struct frame* f3 = apply_transform(Tr2,f2);
      free_transform(Tr2);
      struct transform* Tbp2 = get_random_transform(rna_lib->bp);
      struct transform* Trbp2 = reverse_transform(Tbp2);
      struct frame* f4 = apply_transform(Trbp2,f3);
      struct frame* f12 = average_frames(f1,f2);
      struct frame* f34 = average_frames(f3,f4);
      free_transform(Trbp2);
      set_frame(helix,6*Nbp,f1);
      set_frame(helix,6*Nbp+1,f2);
      set_frame(helix,6*Nbp+2,f12);
      set_frame(helix,6*Nbp+3,f3);
      set_frame(helix,6*Nbp+4,f4);
      set_frame(helix,6*Nbp+5,f34);
      Nbp++;
    }
  }
  return helix;
}

struct trace* build_ideal_forward_helix(struct steptypes* steps, struct trace* trace, struct library* helix_lib) {
  /*
                 Tbp2     
              f3  ->  f4
              ^       |
          Tf2 |  Tbp1 V
              f1  ->  f2
          Tf1 ^       | Tr
              |       V
      fi-1 -> fi  ->  fi+1 -> ...
        T[i-1]   T[i]    T[i+1]
  */
  assert(steps);
  assert(trace);
  assert(helix_lib);
  assert(helix_lib->bb);
  assert(helix_lib->bp);
  int Nbp = steptypes_bp_count(steps);
  struct trace* helix = create_trace(Nbp*6);
  Nbp = 0;
  for(int i=0; i<steptypes_length(steps); i++) {
    char* step = steptypes_get(steps,i);
    if(strcmp(step,"BP") == 0) {
      struct frame* f = get_frame(trace,i);
      struct frame* f_ = get_frame(trace,i+1);
      struct transform* Tf1 = get_transform(helix_lib->bb,0);
      struct frame* f1 = apply_transform(Tf1,f);
      struct transform* Tr1 = reverse_transform(Tf1);
      struct frame* f2 = apply_transform(Tr1,f_);
      struct transform* Tf2 = get_transform(helix_lib->bb,1);
      struct frame* f3 = apply_transform(Tf2,f1);
      struct transform* Tr2 = reverse_transform(Tf2);
      struct frame* f4 = apply_transform(Tr2,f2);
      struct frame* f12 = average_frames(f1,f2);
      struct frame* f34 = average_frames(f3,f4);
      // memory clean-up
      free_transform(Tr1);
      free_transform(Tr2);
      // add the helix
      set_frame(helix,6*Nbp,f1);
      set_frame(helix,6*Nbp+1,f2);
      set_frame(helix,6*Nbp+2,f12); // pseudoatom
      set_frame(helix,6*Nbp+3,f3);
      set_frame(helix,6*Nbp+4,f4);
      set_frame(helix,6*Nbp+5,f34); // pseudoatom
      Nbp++;
    }
  }
  return helix;
}

struct trace* build_ideal_reverse_helix(struct steptypes* steps, struct trace* trace, struct library* helix_lib) {
  /*
                 Tbp2     
              f4  ->  f3
          Tf2 ^       | Tr2
              |  Tbp1 V
              f2  ->  f1
          Tf1 ^       | Tr1
              |       V
      fi-1 -> fi  ->  fi+1 -> ...
        T[i-1]   T[i]    T[i+1]
                  BP
  */
  assert(steps);
  assert(trace);
  assert(helix_lib);
  assert(helix_lib->bb);
  assert(helix_lib->bp);
  int Nbp = steptypes_bp_count(steps);
  struct trace* helix = create_trace(Nbp*6);
  Nbp = 0;
  if(steptypes_length(steps)+1 != trace_length(trace)) {
    printf("error in lengths in build_ideal_reverse_helix: %d vs %d\n",steptypes_length(steps),trace_length(trace));
    exit(1);
  }
  for(int i=steptypes_length(steps)-1; i>=0; i--) {
    char* step = steptypes_get(steps,i);
    if(strcmp(step,"BP") == 0) {
      struct frame* f = get_frame(trace,steptypes_length(steps)-1-i); // there should be 1 more in trace vs steps
      struct frame* f_ = get_frame(trace,steptypes_length(steps)-i);
      struct transform* Tf1 = get_transform(helix_lib->bb,0);
      struct transform* Tr1 = reverse_transform(Tf1);
      struct frame* f1 = apply_transform(Tr1,f);
      struct frame* f2 = apply_transform(Tf1,f_);
      struct transform* Tf2 = get_transform(helix_lib->bb,1);
      struct transform* Tr2 = reverse_transform(Tf2);
      struct frame* f3 = apply_transform(Tr2,f1);
      struct frame* f4 = apply_transform(Tf2,f2);
      struct frame* f12 = average_frames(f1,f2);
      struct frame* f34 = average_frames(f3,f4);
      // memory clean-up
      free_transform(Tr1);
      free_transform(Tr2);
      // add the helix
      set_frame(helix,6*Nbp,f1);
      set_frame(helix,6*Nbp+1,f2);
      set_frame(helix,6*Nbp+2,f12);
      set_frame(helix,6*Nbp+3,f3);
      set_frame(helix,6*Nbp+4,f4);
      set_frame(helix,6*Nbp+5,f34);
      Nbp++;
    }
  }
  return helix;
}

struct trace* build_forward_helix(struct steptypes* fwd_steps, struct trace* fwd_trace, struct library* rna_lib, struct library* helix_lib) {
  // wrapper function depending on whether ideal helix is read in
  if(helix_lib->bp) {
    return build_ideal_forward_helix(fwd_steps,fwd_trace,helix_lib);
  } else {
    return build_random_forward_helix(fwd_steps,fwd_trace,rna_lib);
  }
}

struct trace* build_reverse_helix(struct steptypes* rev_steps, struct trace* rev_trace, struct library* rna_lib, struct library* helix_lib) {
  // wrapper function depending on whether idea helix is read in
  if(helix_lib->bp) {
    return build_ideal_reverse_helix(rev_steps,rev_trace,helix_lib);
  } else {
    return build_random_reverse_helix(rev_steps,rev_trace,rna_lib);
  }
}

struct trace* get_pts_forward(struct steptypes* steps, struct library* rna_lib, int N_samples) {
  assert(steps);
  assert(rna_lib);
  assert(rna_lib->bb);
  assert(rna_lib->bp);
  struct trace* pts_fwd = create_trace(N_samples); // create the library to store the points
  for(int i=0; i<N_samples; i++) {
    struct trace* trace = get_random_forward_trace(steps,rna_lib);
    struct frame* x_fwd = copy_frame(get_frame(trace,-1)); // get the last point in the trace
    add_frame_euler_vector(x_fwd);
    set_frame(pts_fwd,i,x_fwd);
    free_trace(trace);
  }
  return pts_fwd;
}

struct trace* get_pts_reverse(struct steptypes* steps, struct library* rna_lib, int N_samples) {
  assert(steps);
  assert(rna_lib);
  assert(rna_lib->bb);
  assert(rna_lib->bp);
  struct trace* pts_fwd = get_pts_forward(steps,rna_lib,N_samples); // temporary forward trace to be reversed
  struct trace* pts_rev = create_trace(N_samples); // data structure to store reverse trace
  for(int i=0; i<trace_length(pts_fwd); i++) {
    struct frame* x_fwd = get_frame(pts_fwd,i); // T goes from start of steps to end for each replicate
    struct frame* x_rev = reverse_frame(x_fwd); // reverse it to go from end to start
    add_frame_euler_vector(x_rev);
    set_frame(pts_rev,i,x_rev);
  }
  free_trace(pts_fwd); // this should free the individual transforms in pts_fwd
  return pts_rev;
}

struct trace* get_pts_reverse2(struct steptypes* steps, struct library* rna_lib, int N_samples) {
  // alternative approach, not currently used
  assert(steps);
  assert(rna_lib);
  assert(rna_lib->bb);
  assert(rna_lib->bp);
  struct trace* pts_rev = create_trace(N_samples); // create the library to store the points
  for(int i=0; i<N_samples; i++) {
    struct trace* trace = get_random_reverse_trace(steps,rna_lib); // assumes steps fwd, will flip in the fnctn
    struct frame* x_rev = copy_frame(get_frame(trace,-1)); // get the last point in the trace
    add_frame_euler_vector(x_rev);
    set_frame(pts_rev,i,x_rev);
    free_trace(trace);
  }
  return pts_rev;
}

struct trace* read_pts_files(char* pts_file) {
  // read a set of points defined in a file. Used for debugging, not part of current workflow
  int MAX_LINE_LENGTH = 1000;
  FILE* F = fopen(pts_file,"r");
  if(F == NULL) {
    printf("Error opening the file:\"%s\".\n Check permissions/file name.\n",pts_file);
    exit(1);
  }

  struct trace* pts = create_trace(100);
  char line[MAX_LINE_LENGTH];
  int count = 0;
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    struct split_text* row = split_line(line,"\t");
    if(get_split_text_cols(row) == 7) {
      if(strcmp(get_split_text_term(row,0),"idx") != 0) {
	struct frame* f = create_frame(1.0);
	double* ctr = create_vector();
	set_vector_elements(ctr,atof(get_split_text_term(row,1)),atof(get_split_text_term(row,2)),atof(get_split_text_term(row,3)));
	set_frame_center(f,ctr);
	free_vector(ctr);	
	double* EV = create_vector();
	set_vector_elements(EV,atof(get_split_text_term(row,4)),atof(get_split_text_term(row,5)),atof(get_split_text_term(row,6)));
	set_frame_euler_vector(f,EV);
	free_vector(EV);
	set_frame(pts,count,f);
	count++;
      }
    }
  }
  return pts;
}

/*************
 * KDE tools *
 *************/

struct array* get_kde_bandwidth(struct trace* pts) {
  assert(pts);
  // calculate bandwidth using Silverman's rule
  int d = 6; // assuming 6D coordinates
  int n = trace_length(pts); // number of points
  
  struct array* h = create_array(d);
  double* std_dev = get_T6_std_dev(pts);
  
  for(int i=0; i<d; i++) {
    // the h[i] terms correspond to diagonal terms sqrt(H[i,i]) of a bandwidth matrix
    if(std_dev[i] == 0.0) {
      printf("zero std-dev found %d\n",i);
      exit(0);
    }
    array_set(h, i, pow(4.0/(n*(d+2)),1.0/(d+4)) * std_dev[i] );
  }  
  free(std_dev);
  return h;
}

double* get_T6_std_dev_flat(struct array* pts, int N, double* std_dev) {
  // calculate bandwidth using Silverman's rule
  int d = 6; // assuming 6D coordinates
  double sum[d];
  double sum_squares[d];
  
  for(int j=0; j<d; j++) {
    sum[j] = 0.0;
    sum_squares[j] = 0.0;
  }
  for(int i=0; i<N; i++) {
    for(int j=0; j<d; j++) {
      double val = pts->data[i*6 + j];
      sum[j] += val;
      sum_squares[j] += val*val;	
    }
  } 
  for(int j=0; j<d; j++) {
    std_dev[j] = sqrt((sum_squares[j] - sum[j]*sum[j]/N)/(N-1));
  }
}

struct array* get_kde_bandwidth_flat(struct array* pts, int N, double* h) {
  assert(pts);
  // calculate bandwidth using Silverman's rule
  int d = 6; // assuming 6D coordinates
  
  double std_dev[d];
  get_T6_std_dev_flat(pts,N,std_dev);
  //printf("std dev: %f %f %f\n",std_dev[0],std_dev[1],std_dev[3]);
  
  for(int i=0; i<d; i++) {
    // the h[i] terms correspond to diagonal terms sqrt(H[i,i]) of a bandwidth matrix
    if(std_dev[i] == 0.0) {
      printf("zero std-dev found %d\n",i);
      exit(0);
    }
    h[i] = pow(4.0/(N*(d+2)),1.0/(d+4)) * std_dev[i];
  }  
}

double gaussian_kernel(double x) {
  return exp(-0.5*x*x)/sqrt(2.0*pi);
}

double get_circ_prob(struct frame* x_fwd, struct frame* x_rev, struct array* bandwidth) {
  assert(x_fwd);
  assert(x_rev);
  
  int d = 6;
  double product = 1.0;
  double det_bw = 1.0;
  for(int j=0; j<d; j++) {
    double fwd_val = get_T6_value(x_fwd,j);
    double rev_val = get_T6_value(x_rev,j);
    double h = array_get(bandwidth,j);
    double z = (fwd_val - rev_val)/h;
    det_bw *= h;
    product *= gaussian_kernel(z);
  }
  product /= det_bw;
  return product;
}

double get_circ_prob_flat(struct array* pts1, int idx1, struct array* pts2, int idx2, double* bandwidth) {
  assert(pts1);
  assert(pts2);
  assert(bandwidth);
  
  int d=6;
  int offset1 = idx1 * d;
  int offset2 = idx2 * d;
  double product = 1.0;
  double det_bw = 1.0;
  for(int j=0; j<d; j++) {
    double x_val = pts1->data[offset1 + j];
    double f_val = pts2->data[offset2 + j];
    double z = (x_val - f_val)/bandwidth[j];
    product *= gaussian_kernel(z);
    det_bw *= bandwidth[j];
  }
  product /= det_bw;
  return product;
}

double get_kde_probability(struct frame* x, struct trace* pts, struct array* bandwidth) {
  // evaluate the probability f(x) of point x using the data in pts
  // essentially follows https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation
  // using standard multivariate normal kernels
  assert(x);
  assert(pts);
  assert(bandwidth);
  int d = 6;
  int n = trace_length(pts);
  
  // bandwidth is an array, terms of diagonal H matrix
  // determinant is the product of the diagonal terms
  double det_bw = 1.0;
  for(int j=0; j<d; j++) {
    det_bw *= array_get(bandwidth,j);
  }

  double sum = 0.0;
  for(int i=0; i<n; i++) {
    double product = 1.0;
    struct frame* f = get_frame(pts,i);
    for(int j=0; j<d; j++) {
      double x_val = get_T6_value(x,j);
      double f_val = get_T6_value(f,j);
      double h = array_get(bandwidth,j);
      double z = (x_val - f_val)/h;
      product *= gaussian_kernel(z);
    }    
    sum += product;
  }
  sum /= det_bw;
  return sum/(double)n;
}

double get_kde_probability_flat(int idx, struct array* pts1, struct array* pts2, double* bandwidth, int N2) {
  // evaluate the probability f(x) of point x using the data in pts
  // essentially follows https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation
  // using standard multivariate normal kernels
  assert(pts1);
  assert(pts2);
  assert(bandwidth);
  int d = 6;
  
  // bandwidth is an array, terms of diagonal H matrix
  // determinant is the product of the diagonal terms
  double det_bw = 1.0;
  for(int j=0; j<d; j++) {
    det_bw *= bandwidth[j];
  }

  double sum = 0.0;
  int offset1 = idx * d;
  for(int i=0; i<N2; i++) {    
    double product = 1.0;
    int offset2 = i * d;
    for(int j=0; j<d; j++) {
      double x_val = pts1->data[offset1 + j];
      double f_val = pts2->data[offset2 + j];
      double z = (x_val - f_val)/bandwidth[j];
      product *= gaussian_kernel(z);
    }    
    sum += product;
  }
  sum /= det_bw;
  return sum/(double)N2;
}

double get_kde_probability_sap(struct frame* x, struct trace* pts, struct array* bandwidth, struct matrix* M, int k) {
  // sap = self-avoiding polymer: some colliding examples are excluded using M
  // evaluate the probability f(x) of point x using the data in pts
  // essentially follows https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation
  // using standard multivariate normal kernels
  // input int k can either b row or column of input frame, matrix transpose attribute tells which
  assert(x);
  assert(pts);
  assert(bandwidth);
  assert(M);

  int d = 6;
  int n = 0;
  int N_samples = trace_length(pts);
  
  // bandwidth is an array, terms of diagonal H matrix
  // determinant is the product of the diagonal terms
  double det_bw = 1.0;
  for(int j=0; j<d; j++) {
    det_bw *= array_get(bandwidth,j);
  }

  double sum = 0.0;
  for(int i=0; i<N_samples; i++) {
    if(get_element(M,k,i)) {
      double product = 1.0;
      struct frame* f = get_frame(pts,i);
      for(int j=0; j<d; j++) {
	double x_val = get_T6_value(x,j);
	double f_val = get_T6_value(f,j);
	double h = array_get(bandwidth,j);
	double z = (x_val - f_val)/h;
	product *= gaussian_kernel(z);
      }    
      sum += product;
      n++; // should be the total number of self-avoiding pts
    }
  }
  sum /= det_bw;
  if(n == 0) {
    return -1.0;
  } else {
    return sum/(double)n;
  }
}

double get_kde_probability_sap_flat(int idx, struct array* pts1, struct array* pts2, double* bandwidth, struct matrix* M, int N2) {
  // sap = self-avoiding polymer: some colliding examples are excluded using M
  // evaluate the probability f(x) of point x using the data in pts
  // essentially follows https://en.wikipedia.org/wiki/Multivariate_kernel_density_estimation
  // using standard multivariate normal kernels
  // input int k can either b row or column of input frame, matrix transpose attribute tells which
  assert(pts1);
  assert(pts2);
  assert(bandwidth);
  assert(M);
  int d = 6;

  // bandwidth is an array, terms of diagonal H matrix
  // determinant is the product of the diagonal terms
  double det_bw = 1.0;
  for(int j=0; j<d; j++) {
    det_bw *= bandwidth[j];
  }

  int n = 0;
  double sum = 0.0;
  int offset1 = idx * d;
  for(int i=0; i<N2; i++) {
    if(get_element(M,idx,i)) {
      double product = 1.0;
      int offset2 = i * d;
      for(int j=0; j<d; j++) {
	double x_val = pts1->data[offset1 + j];
	double f_val = pts2->data[offset2 + j];
	double z = (x_val - f_val)/bandwidth[j];
	product *= gaussian_kernel(z);
      }
      sum += product;
      n++; // should be the total number of self-avoiding pts
    }
  }
  sum /= det_bw;
  if(n == 0) {
    return -1.0;
  } else {
    return sum/(double)n;
  }
}

double get_rna_closest_approach(struct rna_data* RNA, double bin_size) {  
  // for each RNA chain, compute the closest monomer-monomer distnace
  double d_min = 100000000.0;
  printf("Computing RNA closest approach.\n");
  for(int i=0; i<get_model_count(RNA); i++) {
    struct rna_model* RNA_model = get_rna_model(RNA,i); 
    struct chain* chain = get_first_chain(RNA_model);
    while(chain) {
      int length = chain_length(chain);
      struct trace* chain_trace = create_trace(length);
      struct nuc* n = get_first_nuc(chain);
      int count = 0;
      while(before_end(chain,n)) { // while n is not a sentinel & is a valid nt
	struct frame* f = compute_frame(n);
	// store backbone in trace
	set_frame(chain_trace,count,f);
	count++;
	n = next_nuc(chain,n);
      }
      // at this point, we should examine the chain_trace. first build a detector
      printf("building detector.\n");
      struct detector* detector = create_detector(chain_trace,bin_size);
      printf("allocated memory for detector. %d monomers\n",trace_length(chain_trace));
      for(int i=0; i<trace_length(chain_trace); i++) {
	struct frame* f = get_frame(chain_trace,i);
	double* x = get_frame_center(f);
	int num = get_frame_num(f);
	struct monomer* m = create_monomer(x,num);
	//printf("monomer %d: %.3f %.3f %.3f\n",i,m->r[0],m->r[1],m->r[2]);
	add_monomer_to_detector(detector,m);
      }
      printf("completed detector. Now computing closest approach.\n");
      // print closest approach distances. must be less than bin_size
      char* output_file = "closest_approach.txt";
      FILE* CLOSEST = fopen(output_file,"w");
      for(int i=0; i<trace_length(chain_trace); i++) {
	struct frame* f = get_frame(chain_trace,i);
	double* x = get_frame_center(f);
	int num = get_frame_num(f);
	struct monomer* m = create_monomer(x,num);
	double d = get_closest_approach(detector,m);
	if(d < d_min) {
	  d_min = d;
	}
	if(d < bin_size) {
	  fprintf(CLOSEST,"%.3f\n",d);	  
	}
      }
      fclose(CLOSEST);
      printf("complete. Results stored at %s\n",output_file);
      // move to the next chain if any
      chain = next_chain(chain);
    }
  }
  return d_min;
}

/*********************************
 * Effective Concentration C_eff * 
 *********************************/

double norm_euler_flat(int idx, struct array* pts) {
  int offset = idx*6;
  double theta1 = pts->data[offset + 3];
  double theta2 = pts->data[offset + 4];
  double theta3 = pts->data[offset + 5];
  return sqrt(theta1*theta1 + theta2*theta2 + theta3*theta3);
}

double get_frame_dist_sq_flat(struct array* pts_fwd, int i, struct array* pts_rev, int j) {
  double sum = 0.0;
  for(int k=0; k<3; k++) {
    sum += (pts_fwd->data[i*6 + k] - pts_rev->data[j*6 + k])*(pts_fwd->data[i*6 + k] - pts_rev->data[j*6 + k]);
  }
  return sum;
}

double get_frame_angle_sq_flat(struct array* pts_fwd, int i, struct array* pts_rev, int j) {
  double sum = 0.0;
  for(int k=3; k<6; k++) {
    sum += (pts_fwd->data[i*6 + k] - pts_rev->data[j*6 + k])*(pts_fwd->data[i*6 + k] - pts_rev->data[j*6 + k]);
  }
  return sum;
}

double get_frame_angle_sq_quat(struct array* pts_fwd, int i, struct array* pts_rev, int j) {
  double sum = 0.0;
  double q1[4];
  double q2[4];
  double theta1 = norm_euler_flat(i,pts_fwd);
  if(theta1 < 1e-9) {
    q1[0] = 1.0; q1[1] = 0.0; q1[2] = 0.0; q1[3] = 0.0;
  } else {
    double half_theta = theta1/2.0;
    double s = sin(half_theta);
    double scale = s / theta1;
    q1[0] = cos(half_theta);
    q1[1] = pts_fwd->data[i*6 + 3] * scale;
    q1[2] = pts_fwd->data[i*6 + 4] * scale;
    q1[3] = pts_fwd->data[i*6 + 5] * scale;
  }
  double theta2 = norm_euler_flat(j,pts_rev);
  if(theta2 < 1e-9) {
    q2[0] = 1.0; q2[1] = 0.0; q2[2] = 0.0; q2[3] = 0.0;
  } else {
    double half_theta = theta2/2.0;
    double s = sin(half_theta);
    double scale = s / theta2;
    q2[0] = cos(half_theta);
    q2[1] = pts_rev->data[i*6 + 3] * scale;
    q2[2] = pts_rev->data[i*6 + 4] * scale;
    q2[3] = pts_rev->data[i*6 + 5] * scale;
  }
  double dot = q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2] + q1[3]*q2[3];
  if(dot > 1.0) dot = 1.0;
  if(dot < -1.0) dot = -1.0;
  double angle = 2.0 * acos(fabs(dot));
  return angle*angle;
}


double get_avg_dist_sq(struct array* pts_fwd, struct array* pts_rev, int N) {
  assert(pts_fwd);
  assert(pts_rev);
  
  int sample_count = 0;
  double sum = 0.0;
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      double dist_sq = get_frame_dist_sq_flat(pts_fwd,i,pts_rev,j);
      sum += dist_sq;
      sample_count++;
    }
  }

  // check that data is non-zero
  if(sample_count == 0) {
    printf("Error. No values in C_eff_samples.\n");
    exit(1);
  }

  return sum / sample_count;
}

double get_avg_dist_sq_sap(struct array* pts_fwd, struct array* pts_rev, struct matrix* M, int N1, int N2) {
  assert(pts_fwd);
  assert(pts_rev);
  assert(M);
  
  int sample_count = 0;
  double sum = 0.0;
  for(int i=0; i<N1; i++) {
    for(int j=0; j<N2; j++) {
      if(get_element(M,i,j)) {      
        double dist_sq = get_frame_dist_sq_flat(pts_fwd,i,pts_rev,j);
        sum += dist_sq;
        sample_count++;
      }
    }
  }

  // check that data is non-zero
  if(sample_count == 0) {
    printf("Error. No values in C_eff_samples.\n");
    exit(1);
  }

  return sum / sample_count;
}

double get_avg_angle_sq(struct array* pts_fwd, struct array* pts_rev, int N) {
  assert(pts_fwd);
  assert(pts_rev);  
  
  int sample_count = 0;
  double sum = 0.0;
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      double angle_sq = get_frame_angle_sq_quat(pts_fwd,i,pts_fwd,j);
      sum += angle_sq;
      sample_count++;
    }
  }
  
  // check that data is non-zero
  if(sample_count == 0) {
    printf("Error. No values in C_eff_samples.\n");
    exit(1);
  }
  
  return sum / sample_count;
}

double get_avg_angle_sq_sap(struct array* pts_fwd, struct array* pts_rev, struct matrix* M, int N1, int N2) {
  assert(pts_fwd);
  assert(pts_rev);
  assert(M);
  
  int sample_count = 0;
  double sum = 0.0;
  for(int i=0; i<N1; i++) {
    for(int j=0; j<N2; j++) {
      if(get_element(M,i,j)) {
        double angle_sq = get_frame_angle_sq_quat(pts_fwd,i,pts_fwd,j);
        sum += angle_sq;
        sample_count++;
      }
    }
  }

  // check that data is non-zero
  if(sample_count == 0) {
    printf("Error. No values in C_eff_samples.\n");
    exit(1);
  }
  
  return sum / sample_count;
}

void get_C_eff_samples(struct array* pts_fwd, struct array* pts_rev, double* h_rev, int N_samples, double* C_eff_samples) {
  assert(pts_fwd);
  assert(pts_rev);
  assert(h_rev);
  assert(C_eff_samples);
  
  // pre-calculate multiplier
  const double constant_mult = (8.0 * pi * pi) * (1.0e27 / 6.022e23);
  
  for(int i=0; i<N_samples; i++) {
    double p = get_kde_probability_flat(i, pts_fwd, pts_rev, h_rev, N_samples);
    //printf("%d\t%.3e\n",i,p);
    double v = norm_euler_flat(i, pts_fwd);
    double w = 1.0;
    if(v != 0.0) {
      double sqrt_w = sin(v/2)/(v/2);
      w = sqrt_w*sqrt_w; // more efficient than pow()
    }
    double C_eff = (p / w) * constant_mult;
    C_eff_samples[i] = C_eff;
  }
}

int get_C_eff_samples_sap(struct array* pts_fwd, struct array* pts_rev, double* h_rev, struct matrix* M, int N_fwd, int N_rev, double* C_eff_samples) {
  assert(pts_fwd);
  assert(pts_rev);
  assert(h_rev);
  assert(M);
  assert(C_eff_samples);
  
  double temp_samples[N_fwd];
  // pre-calculate multiplier
  const double constant_mult = (8.0 * pi * pi) * (1.0e27 / 6.022e23);

  int count = 0;
  for(int i=0; i<N_fwd; i++) {
    double p = get_kde_probability_sap_flat(i, pts_fwd, pts_rev, h_rev, M, N_rev);
    if(p != -1) {
      double v = norm_euler_flat(i, pts_fwd);
      double w = 1.0;
      if(v != 0.0) {
        double sqrt_w = sin(v/2)/(v/2);
        w = sqrt_w*sqrt_w; // more efficient than pow()
      }
      double C_eff = (p / w) * constant_mult;
      temp_samples[count] = C_eff;
      count++;
    }
  }
  
  for(int i=0; i<count; i++) {
    C_eff_samples[i] = temp_samples[i];
  }
  return count;
}

struct array* get_C_eff_both_dir(struct steptypes* steps, struct library* rna_lib, int N_samples, double* C_eff_stats) {
  assert(steps);
  assert(rna_lib);
  assert(rna_lib->bb);
  assert(rna_lib->bp);
  int length = steptypes_length(steps);
  int halfway = length/2; // should be floored b/c ints
  int d = 6;

  // forward points
  struct steptypes* fwd_steps = steptypes_splice(steps,0,halfway-1);
  struct array* fwd_pts = create_array(N_samples * d);
  for(int i=0; i<N_samples; i++) {
    struct trace* fwd_trace = get_random_forward_trace(fwd_steps,rna_lib);
    struct frame* x_fwd = get_frame(fwd_trace,-1); // get the last point in the trace
    add_frame_euler_vector(x_fwd);
    for(int j=0; j<6; j++) {
      fwd_pts->data[i*d + j] = get_T6_value(x_fwd,j);
    }
    free_trace(fwd_trace);
  }

  // reverse points
  struct steptypes* rev_steps = steptypes_splice(steps,halfway,length-1);
  struct array* rev_pts = create_array(N_samples * d);
  for(int i=0; i<N_samples; i++) {
    struct trace* rev_trace = get_random_reverse_trace(rev_steps,rna_lib);
    struct frame* x_rev = get_frame(rev_trace,-1); // get the last point in the trace
    add_frame_euler_vector(x_rev);
    for(int j=0; j<d; j++) {
      rev_pts->data[i*d + j] = get_T6_value(x_rev,j);
    }
    free_trace(rev_trace);
  }

  double dist_sq_mean = get_avg_dist_sq(fwd_pts,rev_pts,N_samples);
  double theta_sq_mean = get_avg_angle_sq(fwd_pts,rev_pts,N_samples);
  
  // calculate and store bandwidths for KDE calculations
  double h_fwd[d];
  get_kde_bandwidth_flat(fwd_pts,N_samples,h_fwd);
  double h_rev[d];
  get_kde_bandwidth_flat(rev_pts,N_samples,h_rev);

  double C_eff_fwd_samples[N_samples];
  get_C_eff_samples(fwd_pts,rev_pts,h_rev,N_samples,C_eff_fwd_samples);
  double C_eff_rev_samples[N_samples];
  get_C_eff_samples(rev_pts,fwd_pts,h_fwd,N_samples,C_eff_rev_samples);
  
  double C_eff_fwd_mean = get_sample_mean(C_eff_fwd_samples, N_samples);
  double C_eff_fwd_std_dev = get_sample_std_dev(C_eff_fwd_samples, N_samples);
  double C_eff_fwd_error = C_eff_fwd_std_dev / sqrt(N_samples);
  double C_eff_rev_mean = get_sample_mean(C_eff_rev_samples, N_samples);
  double C_eff_rev_std_dev = get_sample_std_dev(C_eff_rev_samples, N_samples);
  double C_eff_rev_error = C_eff_rev_std_dev / sqrt(N_samples);
  
  // statistical error
  double C_eff_fwd_relerr = C_eff_fwd_error / C_eff_fwd_mean;
  double C_eff_rev_relerr = C_eff_rev_error / C_eff_rev_mean;
  // systematic error
  double C_eff_diff = 0.5 * fabs(log(C_eff_fwd_mean/C_eff_rev_mean));

  // geometric mean of the fwd and rev arithmetic means:
  double C_eff_mean = get_geometric_mean(C_eff_fwd_mean,C_eff_rev_mean);
  // compute relative error as a Euclidian distance of statistical and systematic error "dimensions"
  double C_eff_relerr = sqrt(C_eff_fwd_relerr*C_eff_fwd_relerr + C_eff_rev_relerr*C_eff_rev_relerr + C_eff_diff*C_eff_diff);
  // compute error by scaling unitless relative error by the mean
  double C_eff_error = C_eff_relerr * C_eff_mean;

  // set data to be returned 
  C_eff_stats[0] = C_eff_mean;
  C_eff_stats[1] = C_eff_error;
  C_eff_stats[2] = dist_sq_mean;
  C_eff_stats[3] = theta_sq_mean;

  // memory managment
  free_steptypes(fwd_steps);
  free_steptypes(rev_steps);
  free_array(fwd_pts);
  free_array(rev_pts);
}

struct array* get_C_eff_both_dir_sap(struct motif* motif, struct steptypes* steps, struct library* rna_lib, int N_samples, double bin_size, double d_min, struct library* helix_lib, char* sap_report, double* C_eff_stats) {
  //print_steptypes(steps);
  int length = steptypes_length(steps);
  int halfway = length/2; // should be floored b/c ints
  int d = 6;
  
  // forward traces
  struct steptypes* fwd_steps = steptypes_splice(steps,0,halfway-1);
  struct array* fwd_pts = create_array(N_samples * d);

  // collision detection
  struct trace** fwd_trace = malloc(N_samples * sizeof(struct trace*));
  struct trace** fwd_helix = malloc(N_samples * sizeof(struct trace*));
  struct detector** fwd_detector = malloc(N_samples * sizeof(struct detector*));
  int fwd_collide = 0;
  int N_fwd_pre = 0;
  
  for(int i=0; i<N_samples; i++) {
    int fwd_clear = 1;
    struct trace* trace =  get_random_forward_trace(fwd_steps,rna_lib);
    struct trace* helix = build_forward_helix(fwd_steps,trace,rna_lib,helix_lib);
    struct detector* detector = build_detector_with_helix(trace,helix,bin_size);      
    if(self_collision_helix(detector,trace,helix,d_min)) {
      fwd_clear = 0;
      fwd_collide++;
    }
    if(fwd_clear) {
      fwd_trace[N_fwd_pre] = trace;
      fwd_helix[N_fwd_pre] = helix;
      fwd_detector[N_fwd_pre] = detector;
      struct frame* x_fwd = get_frame(trace,-1); // get the last point in the trace
      add_frame_euler_vector(x_fwd);
      for(int j=0; j<d; j++) {
        fwd_pts->data[N_fwd_pre * d + j] = get_T6_value(x_fwd,j);
      }
      N_fwd_pre++;
    } else {
      free_trace(trace);
      free_trace(helix);
      free_detector(detector);
    }
  }
  
  // reverse traces 
  struct steptypes* rev_steps = steptypes_splice(steps,halfway,length-1);
  struct trace** rev_trace = malloc(N_samples * sizeof(struct trace*));
  struct trace** rev_helix = malloc(N_samples * sizeof(struct trace*));
  struct array* rev_pts = create_array(N_samples * d);
  
  // collision detection
  struct detector** rev_detector = malloc(N_samples * sizeof(struct detector*));
  int rev_collide = 0;
  int N_rev_pre = 0;
  for(int i=0; i<N_samples; i++) {
    struct trace* trace = get_random_reverse_trace(rev_steps,rna_lib);
    struct trace* helix = build_reverse_helix(rev_steps,trace,rna_lib,helix_lib);
    struct detector* detector = build_detector_with_helix(trace,helix,bin_size);
    int rev_clear = 1;
    if(self_collision_helix(detector,trace,helix,d_min)) {
      rev_clear = 0;
      rev_collide++;
    }
    if(rev_clear) {
      // only storing if no collision
      rev_trace[N_rev_pre] = trace;
      rev_helix[N_rev_pre] = helix;
      rev_detector[N_rev_pre] = detector;
      struct frame* x_rev = get_frame(trace,-1); // get the last point in the trace
      add_frame_euler_vector(x_rev);
      for(int j=0; j<d; j++) {
        rev_pts->data[N_rev_pre * d + j] = get_T6_value(x_rev,j);
      }
      N_rev_pre++;
    } else {
      free_trace(trace);
      free_trace(helix);
      free_detector(detector);
    }
  }
  
  // check if the pair avoids each other
  struct matrix* M = initialize_matrix(N_fwd_pre,N_rev_pre);
  int count = 0;
  int total = 0;
  for(int i=0; i<N_fwd_pre; i++) {
    for(int j=0; j<N_rev_pre; j++) {
      total++; // total examined i,j pairs 
      if(trace_collision_helix(fwd_detector[i],rev_trace[j],rev_helix[j],d_min)) {
        // if the rev_trace/helix collides with fwd_trace/helix
        set_element(M,i,j,0);
        count++;
      } 
    }
  }      
  //printf("Completed collision matrix %d\n",omp_get_thread_num());
  
  // print to collisions file
  FILE* SAP_REPORT = fopen(sap_report,"a");
  double p = ((double)count)/total;
  fprint_motif(SAP_REPORT,motif,'\t');
  fprint_steptypes(SAP_REPORT,steps);
  fprintf(SAP_REPORT,"\t%.2e\t%d\t%d\t%d\t%d\t%d\n",p,count,total,fwd_collide,rev_collide,N_samples);
  fclose(SAP_REPORT);

  double dist_sq_mean = get_avg_dist_sq_sap(fwd_pts,rev_pts,M,N_fwd_pre,N_rev_pre);
  double theta_sq_mean = get_avg_angle_sq_sap(fwd_pts,rev_pts,M,N_fwd_pre,N_rev_pre);

  // calculate and store bandwidths for KDE calculations
  // bandwidth calc includes collisions -- could remove them
  double h_fwd[d];
  get_kde_bandwidth_flat(fwd_pts,N_fwd_pre,h_fwd);
  double h_rev[d];
  get_kde_bandwidth_flat(rev_pts,N_rev_pre,h_rev);
  //printf("Completed bandwidth %d\n",omp_get_thread_num());
  
  // compute C_eff samples here
  double C_eff_fwd_samples[N_fwd_pre];
  int N_fwd = get_C_eff_samples_sap(fwd_pts,rev_pts,h_rev,M,N_fwd_pre,N_rev_pre,C_eff_fwd_samples);
  transpose_matrix(M); // transpose matrix to access flipped below
  // for reverse case, all values must be flipped in association:
  double C_eff_rev_samples[N_rev_pre];
  int N_rev = get_C_eff_samples_sap(rev_pts,fwd_pts,h_fwd,M,N_rev_pre,N_fwd_pre,C_eff_rev_samples);
  //printf("samples complete %d\n",omp_get_thread_num());
  
  double C_eff_fwd_mean = get_sample_mean(C_eff_fwd_samples, N_fwd);
  double C_eff_fwd_std_dev = get_sample_std_dev(C_eff_fwd_samples,N_fwd);
  double C_eff_fwd_error = C_eff_fwd_std_dev / sqrt(N_fwd);
  double C_eff_rev_mean = get_sample_mean(C_eff_rev_samples, N_rev);
  double C_eff_rev_std_dev = get_sample_std_dev(C_eff_rev_samples, N_rev);
  double C_eff_rev_error = C_eff_rev_std_dev / sqrt(N_rev);
  
  // statistical error
  double C_eff_fwd_relerr = C_eff_fwd_error / C_eff_fwd_mean;
  double C_eff_rev_relerr = C_eff_rev_error / C_eff_rev_mean;
  // systematic error
  double C_eff_diff = 0.5 * fabs(log(C_eff_fwd_mean/C_eff_rev_mean));
  
  // geometric mean of the fwd and rev arithmetic means:
  double C_eff_mean = get_geometric_mean(C_eff_fwd_mean,C_eff_rev_mean);
  // compute relative error as a Euclidian distance of statistical and systematic error "dimensions"
  double C_eff_relerr = sqrt(C_eff_fwd_relerr*C_eff_fwd_relerr + C_eff_rev_relerr*C_eff_rev_relerr + C_eff_diff*C_eff_diff);
  // compute error by scaling unitless relative error by the mean
  double C_eff_error = C_eff_relerr * C_eff_mean;
  
  // set data to be returned 
  C_eff_stats[0] = C_eff_mean;
  C_eff_stats[1] = C_eff_error;
  C_eff_stats[2] = dist_sq_mean;
  C_eff_stats[3] = theta_sq_mean;
  
  // memory managment
  // fwd detectors and traces:
  for(int i=0; i<N_fwd_pre; i++) {
    free_detector(fwd_detector[i]);
    free_trace(fwd_trace[i]);
    free_trace(fwd_helix[i]);
  }
  free(fwd_detector);
  free(fwd_trace);
  free(fwd_helix);
  // rev detectors and traces:
  for(int j=0; j<N_rev_pre; j++) {
    free_detector(rev_detector[j]);
    free_trace(rev_trace[j]);
    free_trace(rev_helix[j]);
  }
  free(rev_detector);
  free(rev_trace);
  free(rev_helix);
  // free custom collision matrix
  delete_matrix(M);
  //  free the other stuff:
  free_steptypes(fwd_steps);
  free_steptypes(rev_steps);
  free_array(fwd_pts);
  free_array(rev_pts);
}

void fprint_loop_trace(FILE* OUT, int *idx, struct motif* motif, struct steptypes* fwd_steps, struct trace* fwd_trace, struct trace* fwd_helix, struct steptypes* rev_steps, struct trace* rev_trace, struct trace* rev_helix, double p_circ) {
  struct steptypes* steps = steptypes_concat(fwd_steps,rev_steps);
  // labels
  fprintf(OUT,"%d\t",*idx);
  fprint_motif(OUT,motif,'\t');
  fprint_steptypes(OUT,steps);
  // forward trace
  fprintf(OUT,"\t%.2e\t[",p_circ);
  for(int i=0; i<steptypes_length(fwd_steps); i++) {
    char* step = steptypes_get(fwd_steps,i);
    struct frame* x1 = get_frame(fwd_trace,i);
    double* x1_ctr = get_frame_center(x1);
    struct frame* x2 = get_frame(fwd_trace,i+1);
    double* x2_ctr = get_frame_center(x2);
    fprintf(OUT,"((%.2f,%.2f,%.2f),(%.2f,%.2f,%.2f),\"%s\")",x1_ctr[0],x1_ctr[1],x1_ctr[2],x2_ctr[0],x2_ctr[1],x2_ctr[2],step);
    if(i < steptypes_length(fwd_steps)-1) {
      fprintf(OUT,",");
    }    
  }
  fprintf(OUT,"],");

  // reverse trace
  fprintf(OUT,"[");
  for(int i=steptypes_length(rev_steps)-1; i>=0; i--) {
    char* step = steptypes_get(rev_steps,i);
    struct frame* x1 = get_frame(rev_trace,steptypes_length(rev_steps)-1-i);
    double* x1_ctr = get_frame_center(x1);
    struct frame* x2 = get_frame(rev_trace,steptypes_length(rev_steps)-i);
    double* x2_ctr = get_frame_center(x2);
    fprintf(OUT,"((%.2f,%.2f,%.2f),(%.2f,%.2f,%.2f),\"%s\")",x1_ctr[0],x1_ctr[1],x1_ctr[2],x2_ctr[0],x2_ctr[1],x2_ctr[2],step);
    if(i > 0) {
      fprintf(OUT,",");
    }    
  }
  fprintf(OUT,"],");

  // fwd helix
  fprintf(OUT,"[");
  for(int i=0; i<trace_length(fwd_helix)-2; i+=3) {
    struct frame* x1 = get_frame(fwd_helix,i);
    double* x1_ctr = get_frame_center(x1);
    struct frame* x2 = get_frame(fwd_helix,i+1);
    double* x2_ctr = get_frame_center(x2);
    struct frame* x3 = get_frame(fwd_helix,i+2);
    double* x3_ctr = get_frame_center(x3);
    fprintf(OUT,"((%.2f,%.2f,%.2f),(%.2f,%.2f,%.2f),\"%s\"),",x1_ctr[0],x1_ctr[1],x1_ctr[2],x2_ctr[0],x2_ctr[1],x2_ctr[2],"BP");
    fprintf(OUT,"((%.2f,%.2f,%.2f),\"%s\")",x3_ctr[0],x3_ctr[1],x3_ctr[2],"PA");
    if(i<trace_length(fwd_helix)-3) {
      fprintf(OUT,",");
    }
  }
  fprintf(OUT,"],");

  // rev helix
  fprintf(OUT,"[");
  for(int i=0; i<trace_length(rev_helix); i+=3) {
    struct frame* x1 = get_frame(rev_helix,i);
    double* x1_ctr = get_frame_center(x1);
    struct frame* x2 = get_frame(rev_helix,i+1);
    double* x2_ctr = get_frame_center(x2);
    struct frame* x3 = get_frame(rev_helix,i+2);
    double* x3_ctr = get_frame_center(x3);
    fprintf(OUT,"((%.2f,%.2f,%.2f),(%.2f,%.2f,%.2f),\"%s\"),",x1_ctr[0],x1_ctr[1],x1_ctr[2],x2_ctr[0],x2_ctr[1],x2_ctr[2],"BP");
    fprintf(OUT,"((%.2f,%.2f,%.2f),\"%s\")",x3_ctr[0],x3_ctr[1],x3_ctr[2],"PA");
    if(i<trace_length(rev_helix)-3) {
      fprintf(OUT,",");
    }
  }
  fprintf(OUT,"]\n");
}

void print_loop_traces(FILE* OUT, int *idx, struct motif* motif, struct steptypes* steps, struct library* rna_lib, int N_samples, struct library* helix_lib, double p_circ_min) {
  //print_steptypes(steps);
  int length = steptypes_length(steps);
  int halfway = length/2; // should be floored b/c ints
  int d = 6;
  
  // forward traces
  struct steptypes* fwd_steps = steptypes_splice(steps,0,halfway-1);
  struct array* fwd_pts = create_array(N_samples * d);
  struct trace** fwd_trace = malloc(N_samples * sizeof(struct trace*));
  struct trace** fwd_helix = malloc(N_samples * sizeof(struct trace*));
  
  // forward points
  for(int i=0; i<N_samples; i++) {
    int fwd_clear = 1;
    struct trace* trace =  get_random_forward_trace(fwd_steps,rna_lib);
    struct trace* helix = build_forward_helix(fwd_steps,trace,rna_lib,helix_lib);
    struct frame* x_fwd = get_frame(trace,-1); // get the last point in the trace
    add_frame_euler_vector(x_fwd);
    for(int j=0; j<d; j++) {
      fwd_pts->data[i * d + j] = get_T6_value(x_fwd,j);
    }
    fwd_trace[i] = trace;
    fwd_helix[i] = helix;
  }
  
  // reverse traces 
  struct steptypes* rev_steps = steptypes_splice(steps,halfway,length-1);
  struct trace** rev_trace = malloc(N_samples * sizeof(struct trace*));
  struct trace** rev_helix = malloc(N_samples * sizeof(struct trace*));
  struct array* rev_pts = create_array(N_samples * d);

  // reverse points 
  for(int i=0; i<N_samples; i++) {
    struct trace* trace = get_random_reverse_trace(rev_steps,rna_lib);
    struct trace* helix = build_reverse_helix(rev_steps,trace,rna_lib,helix_lib);
    struct frame* x_rev = get_frame(trace,-1); // get the last point in the trace
    add_frame_euler_vector(x_rev);
    for(int j=0; j<d; j++) {
      rev_pts->data[i * d + j] = get_T6_value(x_rev,j);
    }
    rev_trace[i] = trace;
    rev_helix[i] = helix;
  }
  
  // calculate and store bandwidths for KDE calculations
  // bandwidth calc includes collisions -- could remove them
  double h_fwd[d];
  get_kde_bandwidth_flat(fwd_pts,N_samples,h_fwd);
  double h_rev[d];
  get_kde_bandwidth_flat(rev_pts,N_samples,h_rev);
  
  for(int idx1=0; idx1<N_samples; idx1++) {
    for(int idx2=0; idx2<N_samples; idx2++) {
      double p_circ = get_circ_prob_flat(fwd_pts,idx1,rev_pts,idx2,h_rev);
      if(p_circ >= p_circ_min) {
        fprint_loop_trace(OUT,idx,motif,fwd_steps,fwd_trace[idx1],fwd_helix[idx1],rev_steps,rev_trace[idx2],rev_helix[idx2],p_circ);
        (*idx)++;
      }
    }
  }
 
  // memory managment
  for(int i=0; i<N_samples; i++) {
    free_trace(fwd_trace[i]);
    free_trace(fwd_helix[i]);
  }
  free(fwd_trace);
  free(fwd_helix);
  // rev detectors and traces:
  for(int j=0; j<N_samples; j++) {
    free_trace(rev_trace[j]);
    free_trace(rev_helix[j]);
  }
  free(rev_trace);
  free(rev_helix);
  //  free the other stuff:
  free_steptypes(fwd_steps);
  free_steptypes(rev_steps);
  free_array(fwd_pts);
  free_array(rev_pts);
}

void print_loop_traces_sap(FILE* OUT, int *idx, struct motif* motif, struct steptypes* steps, struct library* rna_lib, int N_samples, double bin_size, double d_min, struct library* helix_lib, char* sap_report, double p_circ_min) {
  //print_steptypes(steps);
  int length = steptypes_length(steps);
  int halfway = length/2; // should be floored b/c ints
  int d = 6;
  
  // forward traces
  struct steptypes* fwd_steps = steptypes_splice(steps,0,halfway-1);
  struct array* fwd_pts = create_array(N_samples * d);

  // collision detection
  struct trace** fwd_trace = malloc(N_samples * sizeof(struct trace*));
  struct trace** fwd_helix = malloc(N_samples * sizeof(struct trace*));
  struct detector** fwd_detector = malloc(N_samples * sizeof(struct detector*));
  int fwd_collide = 0;
  int N_fwd_pre = 0;
  
  for(int i=0; i<N_samples; i++) {
    int fwd_clear = 1;
    struct trace* trace =  get_random_forward_trace(fwd_steps,rna_lib);
    struct trace* helix = build_forward_helix(fwd_steps,trace,rna_lib,helix_lib);
    struct detector* detector = build_detector_with_helix(trace,helix,bin_size);      
    if(self_collision_helix(detector,trace,helix,d_min)) {
      fwd_clear = 0;
      fwd_collide++;
    }
    if(fwd_clear) {
      fwd_trace[N_fwd_pre] = trace;
      fwd_helix[N_fwd_pre] = helix;
      fwd_detector[N_fwd_pre] = detector;
      struct frame* x_fwd = get_frame(trace,-1); // get the last point in the trace
      add_frame_euler_vector(x_fwd);
      for(int j=0; j<d; j++) {
        fwd_pts->data[N_fwd_pre * d + j] = get_T6_value(x_fwd,j);
      }
      N_fwd_pre++;
    } else {
      free_trace(trace);
      free_trace(helix);
      free_detector(detector);
    }
  }
  
  // reverse traces 
  struct steptypes* rev_steps = steptypes_splice(steps,halfway,length-1);
  struct trace** rev_trace = malloc(N_samples * sizeof(struct trace*));
  struct trace** rev_helix = malloc(N_samples * sizeof(struct trace*));
  struct array* rev_pts = create_array(N_samples * d);
  
  // collision detection
  struct detector** rev_detector = malloc(N_samples * sizeof(struct detector*));
  int rev_collide = 0;
  int N_rev_pre = 0;
  for(int i=0; i<N_samples; i++) {
    struct trace* trace = get_random_reverse_trace(rev_steps,rna_lib);
    struct trace* helix = build_reverse_helix(rev_steps,trace,rna_lib,helix_lib);
    struct detector* detector = build_detector_with_helix(trace,helix,bin_size);
    int rev_clear = 1;
    if(self_collision_helix(detector,trace,helix,d_min)) {
      rev_clear = 0;
      rev_collide++;
    }
    if(rev_clear) {
      // only storing if no collision
      rev_trace[N_rev_pre] = trace;
      rev_helix[N_rev_pre] = helix;
      rev_detector[N_rev_pre] = detector;
      struct frame* x_rev = get_frame(trace,-1); // get the last point in the trace
      add_frame_euler_vector(x_rev);
      for(int j=0; j<d; j++) {
        rev_pts->data[N_rev_pre * d + j] = get_T6_value(x_rev,j);
      }
      N_rev_pre++;
    } else {
      free_trace(trace);
      free_trace(helix);
      free_detector(detector);
    }
  }
  
  // check if the pair avoids each other
  struct matrix* M = initialize_matrix(N_fwd_pre,N_rev_pre);
  for(int i=0; i<N_fwd_pre; i++) {
    for(int j=0; j<N_rev_pre; j++) {
      if(trace_collision_helix(fwd_detector[i],rev_trace[j],rev_helix[j],d_min)) {
        // if the rev_trace/helix collides with fwd_trace/helix
        set_element(M,i,j,0);
      }
    }
  }

  // calculate and store bandwidths for KDE calculations
  // bandwidth calc includes collisions -- could remove them
  double h_fwd[d];
  get_kde_bandwidth_flat(fwd_pts,N_fwd_pre,h_fwd);
  double h_rev[d];
  get_kde_bandwidth_flat(rev_pts,N_fwd_pre,h_rev);

  for(int id1=0; id1<N_fwd_pre; id1++) {
    for(int id2=0; id2<N_rev_pre; id2++) {
      if(get_element(M,id1,id2)) {
        double p_circ = get_circ_prob_flat(fwd_pts,id1,rev_pts,id2,h_rev);
        if(p_circ >= p_circ_min) {
          fprint_loop_trace(OUT,idx,motif,fwd_steps,fwd_trace[id1],fwd_helix[id1],rev_steps,rev_trace[id2],rev_helix[id2],p_circ);
          (*idx)++;
        }
      }
    }
  }
  
  // memory managment
  // fwd detectors and traces:
  for(int i=0; i<N_fwd_pre; i++) {
    free_detector(fwd_detector[i]);
    free_trace(fwd_trace[i]);
    free_trace(fwd_helix[i]);
  }
  free(fwd_detector);
  free(fwd_trace);
  free(fwd_helix);
  // rev detectors and traces:
  for(int j=0; j<N_rev_pre; j++) {
    free_detector(rev_detector[j]);
    free_trace(rev_trace[j]);
    free_trace(rev_helix[j]);
  }
  free(rev_detector);
  free(rev_trace);
  free(rev_helix);
  // free custom collision matrix
  delete_matrix(M);
  //  free the other stuff:
  free_steptypes(fwd_steps);
  free_steptypes(rev_steps);
  free_array(fwd_pts);
  free_array(rev_pts);
}

