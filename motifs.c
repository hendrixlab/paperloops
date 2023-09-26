 /*
 * This library contains the definitions of structures and functions related to loop motifs
 * defining the lengths of each loop in a multiloop structure. e.g. [3,4,2,1]
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#include "motifs.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * data structures for the transform and for the library of transforms
 */

struct motif {
  int* loop;
  int num;
  int id;
  struct motif* next;
};

struct motif_group {
  struct motif** data;
  int count; // current number of motifs
  int size; // should have a defined size
  int N; // number of loops in motif group
  int L; // total length of all loops in each motif
};

struct motif_library {
  struct motif_group** data;
  int N_loops;
  int L_max;
};

/*
 * Tools for building/adding/removing to a motif
 */

struct motif* create_motif(int num) {
  struct motif* motif = malloc(sizeof(struct motif));
  motif->num = num;
  motif->loop = malloc(sizeof(int)*num);
  motif->next = NULL;
  return motif;
}

void initialize_motif(struct motif* motif, int loop) {
  set_motif_loop(motif,0,loop);
  set_motif_id(motif,loop);
}

void set_motif_loop(struct motif* motif, int idx, int loop) {
  assert(motif);
  motif->loop[idx] = loop;
}

void set_motif_id(struct motif* motif, int id) {
  assert(motif);
  motif->id = id;
}

int get_motif_num(struct motif* motif) {
  return motif->num;
}

int get_motif_loop(struct motif* motif, int idx) {
  assert(motif);
  return motif->loop[idx];
}

int get_motif_id(struct motif* motif) {
  return motif->id;
}

int get_motif_sum(struct motif* motif) {
  int sum = 0;
  int num = get_motif_num(motif);
  for(int i=0; i<num; i++) {
    sum += get_motif_loop(motif,i);
  }
  return sum;
}

int print_motif(struct motif* motif, char c) {
  int n = get_motif_num(motif);
  printf("motif %d: (",motif->id);
  for(int i=0; i<n; i++) {
    printf("%d",motif->loop[i]);
    if(i==n-1) {
      printf(")%c",c);
    } else {
      printf(", ");
    }
  }  
}

int get_new_id(int id, int new_num, int new_loop, int L_max) {
  int new_id = id + new_loop*pow(L_max+1,new_num-1); // next place in base-(L_max+1) number  
}

int compute_motif_id(struct motif* motif, int L_max) {
  int new_id = 0;
  int num = get_motif_num(motif);
  for(int i=0; i<num; i++) {
    new_id += get_motif_loop(motif,i)*pow(L_max+1,i);
  }
  return new_id;
}

struct motif* expand_motif(struct motif* motif, int new_loop, int L_max) {
  // adds to a motif a new loop length
  int num = get_motif_num(motif);
  int id = get_motif_id(motif);

  int new_num = num+1;
  struct motif* new_motif = create_motif(new_num);
  int new_id = get_new_id(id,new_num,new_loop,L_max);
  set_motif_id(new_motif,new_id);
  
  for(int i=0;i<num;i++) { // copy the old motif
    int loop = get_motif_loop(motif,i);
    set_motif_loop(new_motif,i,loop);
  }
  set_motif_loop(new_motif,num,new_loop); // add the new loop
  
  return new_motif;
}

struct motif* rotate_motif(struct motif* motif, int L_max) {
  // adds to a motif a new loop length
  int num = get_motif_num(motif);

  struct motif* rot_motif = create_motif(num);
    
  for(int i=1;i<num;i++) { // copy the old motif
    int loop = get_motif_loop(motif,i);
    set_motif_loop(rot_motif,i-1,loop);
  }
  int loop = get_motif_loop(motif,0);
  set_motif_loop(rot_motif,num-1,loop);
  
  int new_id = compute_motif_id(rot_motif,L_max);
  set_motif_id(rot_motif,new_id);

  // important: building in linked list of rotations attached to main motif
  //motif->next = rot_motif;
  
  return rot_motif;
}

struct motif* get_next_motif(struct motif* motif) {
  assert(motif);
  return motif->next;
}

struct motif* set_next_motif(struct motif* motif, struct motif* next) {
  assert(motif);
  return motif->next = next;
}



/*
 * Tools for building/adding/removing elements from motif_group
 */

struct motif_group* create_motif_group(int size) {
  struct motif_group* motif_group = malloc(sizeof(struct motif_group));  
  motif_group->count = 0;  
  motif_group->size = size;
  motif_group->data = malloc(sizeof(struct motif*) * size);
  return motif_group;
}

void free_motif_group(struct motif_group* motif_group) {
  assert(motif_group);
  free(motif_group->data);
  free(motif_group);
}

int get_motif_group_size(struct motif_group* motif_group) {
  assert(motif_group);
  return motif_group->size;
}

int get_motif_group_count(struct motif_group* motif_group) {
  assert(motif_group);
  return motif_group->count;
}

struct motif* get_motif(struct motif_group* motif_group, int idx) {
  assert(motif_group);
  if((idx < motif_group->size && idx >= 0)) {
    return motif_group->data[idx];
  }
  return NULL;
}

int motif_empty(struct motif_group* motif_group, int idx) {
  assert(idx < motif_group->size && idx >= 0);
  if(motif_group->data[idx]) {
    return 0;
  }
  return 1;
}

void add_motif(struct motif_group* motif_group, struct motif* motif) {
  assert(motif_group);
  assert(motif);

  int idx = get_motif_id(motif);
  assert(idx < motif_group->size && idx >= 0);
  
  motif_group->data[idx] = motif;
  motif_group->count++;
}

void append_motif(struct motif_group* motif_group, struct motif* motif) {
  assert(motif_group);
  assert(motif);

  int len = motif_group->count;
  if(len == motif_group->size) {
    // resize the motif_group
    int new_size = 2 * motif_group->size;
    struct motif** new_data = malloc(new_size * sizeof(struct motif*));
    for(int i=0; i<motif_group->count; i++) {
      new_data[i] = motif_group->data[i];
    }
    free(motif_group->data);
    motif_group->data = new_data;
    motif_group->size = new_size;
    //printf("resizing to %d\n",new_size);
  } 
  motif_group->data[len] = motif;
  motif_group->count++;
}

struct motif_group*** create_motif_lib(int N_loops, int L_max) {
  struct motif_group*** motif_lib = malloc((N_loops+1)*sizeof(struct motif_group**));
  for(int n=0;n<=N_loops;n++) {
    motif_lib[n] = malloc((L_max)*sizeof(struct motif_group*));
    for(int L=0; L<=L_max; L++) {
      motif_lib[n][L] = create_motif_group(1);
    }
  }
  return motif_lib;
}

struct motif_group** create_motif_list(int N_loops, int L_max) {
  struct motif_group** motif_list = malloc((N_loops+1)*sizeof(struct motif_group*));  
  motif_list[0] = create_motif_group(0); // sort of a place holder, has size=0
  motif_list[1] = create_motif_group(L_max+1); // hairpin loops, should be L_max in size
  
  for(int L=0;L<=L_max;L++) {
    struct motif* motif = create_motif(1);
    initialize_motif(motif,L); // one-loop motifs, setting l_0 = L;
    add_motif(motif_list[1],motif);
  }
  return motif_list;
}

struct motif_group*** build_motif_list(int N_loops, int L_max) {
  // build motif library
  struct motif_group*** motif_lib = create_motif_lib(N_loops,L_max);
  // build motif_group structure, initialize with singleton motifs:
  struct motif_group** motif_list = create_motif_list(N_loops,L_max);  
  // recursively build motif_lib for different lengths  
  
  for(int n=2;n<=N_loops;n++) { // create multi-junction loops 
    int new_size = 1 + pow(L_max+1,n); // 1 + largest possible id for n loops that add up to L_max or less
    motif_list[n] = create_motif_group(new_size);
    int size = get_motif_group_size(motif_list[n-1]); // size of n-1
    for(int id=0;id<size;id++) { // id of n-1 group
      if(!motif_empty(motif_list[n-1],id)) { // checking if valid id for the n-1 group 
	struct motif* motif = get_motif(motif_list[n-1],id);
	int sum = get_motif_sum(motif);
	for(int L=0;L<=L_max-sum;L++) {
	  int new_id = get_new_id(id,n,L,L_max);
	  if(motif_empty(motif_list[n],new_id)) { // only add the new motif if it is not filled yet.
	    int L_tot = sum + L;
	    struct motif* new_motif = expand_motif(motif,L,L_max);
	    printf("adding motif by expansion: ");
	    print_motif(new_motif,'\n');
	    add_motif(motif_list[n],new_motif);	    
	    append_motif(motif_lib[n][L_tot],new_motif);
	    for(int j=0; j<n; j++) { // apart from repeats, motif with n loops has n rotations
	      struct motif* prev_motif = new_motif;
	      new_motif = rotate_motif(new_motif,L_max);
	      int new_id = get_motif_id(new_motif);
	      if(motif_empty(motif_list[n],new_id)) { // only add if it is not visited/stored yet
		printf("adding motif by rotation: ");
		print_motif(new_motif,'\n');
		add_motif(motif_list[n],new_motif);
		set_next_motif(prev_motif,new_motif);
	      }
 	    }
	  }
	}
      }
    }
  }
  return motif_lib;
}
