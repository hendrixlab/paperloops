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
  if(argc == 4) {
    printf("input PDB file name: %s\n",argv[1]);
    printf("input base pair file name: %s\n",argv[2]);
    printf("input output suffix: %s\n",argv[3]);
  } else {
    printf("Usage: %s <PDB file name> <base pair file>\n",argv[0]);
    exit(1);
  }
  char* pdb_file_name = argv[1];
  char* basepair_file_name = argv[2];
  char* suffix = argv[3];
  
  srand(time(NULL)); // set random seed by time
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
  
  int N_loops = 6;
  int L_max = 20;
  int N_iter = 2000;

  // actual motifs built with efficient algorithm
  struct motif_group*** M = build_motif_lib(N_loops,L_max);
  // theoretical counts based on equation
  int** N = necklace_counts_fixed_length(N_loops,L_max);
 
  char out_all_file[50] = "C_eff_all_results_";
  char out_file[50] = "C_eff_results_";
  strcat(out_all_file,suffix);
  strcat(out_file,suffix);
  strcat(out_all_file,".txt");
  strcat(out_file,".txt");
    
  FILE* OUT_ALL = fopen(out_all_file,"w");
  FILE* OUT = fopen(out_file,"w");
  for(int n=1; n<=N_loops; n++) {
    if(n > 1) {
      L_max = 4;
    }
    for(int L=0; L<=L_max; L++) {
      int count = get_motif_group_count(M[n][L]);
      printf("N_pred(n=%d,L=%d)=%d, N_obs=%d\n",n,L,N[n][L],count);
      fflush(stdout);
      for(int i=0; i<count; i++) {
	struct motif* m = get_motif(M[n][L],i);
	while(m) {	
	  // here we have a motif, part of a necklace, first one should be lowest lexicographically.
	  struct steptypes* steps = motif_to_steptypes(m);	  	  
	  fprint_motif(OUT,m,'\t');
	  //print_motif(m,'\n');
	  //fprint_steptypes(OUT,steps);	  
	  double num = 0.0;
	  double den = 0.0;
	  for(int j=0; j<steptypes_length(steps); j++) {
	    struct array* C_eff_stats = get_C_eff_both_dir(steps,bb_lib,bp_lib,N_iter);
	    double C_eff_mean = array_get(C_eff_stats,0);
	    double C_eff_error = array_get(C_eff_stats,1);
	    fprint_motif(OUT_ALL,m,'\t');
	    fprintf(OUT_ALL,"%.8f\t%.8f\n",C_eff_mean,C_eff_error);
	    double C_eff_relerr = C_eff_error / C_eff_mean;
	    num += log(C_eff_mean)/(C_eff_relerr*C_eff_relerr);
	    den += 1.0/(C_eff_relerr*C_eff_relerr);
	    // rotate the steptypes and add to running average
	    steps = rotate_steptypes(steps);
	    // memory management:
	    free_array(C_eff_stats);
	  }
	  // weighted average over steptypes rotations:
	  double full_C_eff_mean = exp(num/den);
	  double full_C_eff_err = sqrt(1.0/den)*full_C_eff_mean;
	  fprintf(OUT,"%.8f\t%.8f\n",full_C_eff_mean,full_C_eff_err);
	  fflush(OUT);
	  m = get_next_motif(m); // rotate the motif to get an equivalent one in different order
	  free_steptypes(steps);
	}
      }
    }
  }
  free_motif_lib(M,N_loops,L_max);
  free_necklace_counts(N,N_loops);
  free_library(bb_lib);
  free_library(bp_lib);
  return 0;
}
