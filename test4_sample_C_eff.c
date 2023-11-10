#include "rna_io.h"
#include "rna_data.h"
#include "simulation.h"
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
  struct library* bb_lib = build_backbone_transform_library(RNA);
  printf("basepair libraries:\n");
  struct library* bp_lib = build_basepair_transform_library(RNA,basepairs);
  printf("completed.\n");

  int N_loops = 5;
  int L_max = 8;
  int N_iter = 200;
  // actual motifs built with efficient algorithm
  struct motif_group*** M = build_motif_list(N_loops,L_max);
  // theoretical lengths based on equation
  int** N = necklace_counts_fixed_length(N_loops,L_max);
 
  FILE* OUT = fopen("C_eff_results4.txt","w");
  for(int n=1; n<=N_loops; n++) {
    for(int L=0; L<=L_max; L++) {
      int count = get_motif_group_count(M[n][L]);
      printf("N(n=%d,L=%d)=%d, N_obs=%d\n",n,L,N[n][L],count);
      for(int i=0; i<count; i++) {
	// here we have a motif, part of a necklace, first one should be lowest lexicographically.
	struct motif* m = get_motif(M[n][L],i);
	while(m) {	
	  struct steptypes* steps = motif_to_steptypes(m);	  
	  fprint_motif(OUT,m,'\t');
	  fprint_steptypes(OUT,steps);
	  struct array* C_eff_samples = get_C_eff_both_dir(steps,bb_lib,bp_lib,N_iter);
	  double C_eff_mean = get_sample_mean(C_eff_samples);
	  double C_eff_std_dev = get_sample_std_dev(C_eff_samples);
	  fprintf(OUT,"\t%.6f\t%.6f\n",C_eff_mean,C_eff_std_dev);
	  m = get_next_motif(m);
	}
      }
    }
  }
  fclose(OUT);
  return 0;
}
