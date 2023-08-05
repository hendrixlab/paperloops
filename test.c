#include "rna_pdb.h"
#include "rna_data.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  if(argc == 2) {
    printf("input file name: %s\n",argv[1]);
  } else {
    printf("Usage: %s <PDB file name>\n",argv[0]);
    exit(1);
  }
  char* input_file_name = argv[1];

  struct rna_data* RNA = read_pdb_file(input_file_name);
  get_dinucleotides(RNA);
}
