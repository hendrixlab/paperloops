#include "rna_io.h"
#include "rna_data.h"
#include "simulation.h"
#include "transform.h"
#include "motifs.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[]) {
  if(argc == 3) {
    printf("input PDB file name: %s\n",argv[1]);
    printf("input base pair file name: %s\n",argv[2]);
  } else {
    printf("Usage: %s <PDB file name> <base pair file>\n",argv[0]);
    exit(1);
  }
  char* pdb_file_name = argv[1];
  char* basepair_file_name = argv[2];
  
  srand(time(NULL)); // set random seed by time
  //srand(777);
  int r = rand() % 7; // discard first value based on stackoverflow
  
  struct rna_data* RNA = read_pdb_file(pdb_file_name);
  struct bp_array* basepairs = read_basepair_file(basepair_file_name);
  
  printf("Building transform libraries.\n");
  printf("backbone libraries:\n");
  struct library* bb_lib = build_backbone_transform_library(RNA);
  printf("basepair libraries:\n");
  struct library* bp_lib = build_basepair_transform_library(RNA,basepairs);
  printf("libraries complete.\n");
  free_bp_array(basepairs);
  
  char bb_transform_file[50] = "bb_transform_file.txt";
  FILE* BB_OUT = fopen(bb_transform_file,"w");
  char bp_transform_file[50] = "bp_transform_file.txt";
  FILE* BP_OUT = fopen(bp_transform_file,"w");
  
  int d = 6;

  for(int i=0; i<library_size(bb_lib); i++) {
    struct transform* T = get_transform(bb_lib,i);
    double* EV = compute_euler_vector(get_rotation_matrix(T));
    set_euler_vector(T,EV);
    free_vector(EV);
    fprintf(BB_OUT,"%d",i);
    for(int j=0; j<d; j++) {
      double x = get_T6_value(T,j);
      fprintf(BB_OUT,"\t%.6f",x);
    }
    fprintf(BB_OUT,"\n");
  }

  for(int i=0; i<library_size(bp_lib); i++) {
    struct transform* T = get_transform(bp_lib,i);
    double* EV = compute_euler_vector(get_rotation_matrix(T));
    set_euler_vector(T,EV);
    free_vector(EV);
    fprintf(BP_OUT,"%d",i);
    for(int j=0; j<d; j++) {
      double x = get_T6_value(T,j);
      fprintf(BP_OUT,"\t%.6f",x);
    }
    fprintf(BP_OUT,"\n");
  }
  
  free_library(bb_lib);
  free_library(bp_lib);
  return 0;
}
