#include "rna_pdb.h"
#include "rna_data.h"
#include "transform.h"
#include "motifs.h"

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

  int N_loops = 8;
  int L_max = 10;
  // actual motifs built with efficient algorithm
  struct motif_group*** M = build_motif_list(N_loops,L_max);
  // theoretical lengths based on equation
  int** N = necklace_counts_fixed_length(N_loops,L_max);

  for(int n=1; n<=N_loops; n++) {
    for(int L=0; L<=L_max; L++) {
      int count = get_motif_group_count(M[n][L]);
      printf("N(n=%d,L=%d)=%d, count=%d\n",n,L,N[n][L],count);
      for(int i=0; i<count; i++) {
	printf("%d-",i);
	struct motif* m = get_motif(M[n][L],i);
	while(m) {
	  print_motif(m,' ');
	  m = get_next_motif(m);
	}
	printf("\n");
      }
    }
  }
  return 0;
}
