#include "rna_io.h"
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
  printf("library complete.\n");

  int N_loops = 4;
  int L_max = 10;
  // actual motifs built with efficient algorithm
  struct motif_group*** M = build_motif_lib(N_loops,L_max);
  // theoretical lengths based on equation
  int** N = necklace_counts_fixed_length(N_loops,L_max);

  printf("printing out motifs.\n");
  FILE* OUT = fopen("loop_motif_list.txt","w");
  for(int n=1; n<=N_loops; n++) {
    for(int L=0; L<=L_max; L++) {
      int count = get_motif_group_count(M[n][L]);
      fprintf(OUT,"# Npred(n=%d,L=%d)=%d, Nobs=%d\n",n,L,N[n][L],count);
      for(int i=0; i<count; i++) {
	fprintf(OUT,"%d ",i+1);
	struct motif* m = get_motif(M[n][L],i);
	while(m) {
	  fprint_motif(OUT,m,' ');
	  m = get_next_motif(m);
	}
	fprintf(OUT,"\n");
      }
    }
  }
  printf("complete.\n");
  return 0;
}
