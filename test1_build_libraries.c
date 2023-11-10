#include "rna_io.h"
#include "rna_data.h"
#include "transform.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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

  struct rna_data* RNA = read_pdb_file(pdb_file_name);
  struct bp_array* basepairs = read_basepair_file(basepair_file_name);
  
  printf("Building transform libraries.\n");
  printf("backbone libraries:\n");
  build_backbone_transform_library(RNA);
  printf("basepair libraries:\n");
  build_basepair_transform_library(RNA,basepairs);
  printf("completed.\n");

  return 0;
}
