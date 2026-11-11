#include "rna_io.h"
#include "rna_data.h"
#include "simulation.h"
#include "transform.h"
#include "motifs.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>


// struct to store the output buffer
struct output_buffer {
  char** data;
  int len;
  int count;
};

/**********************
 * FUNCTION TEMPLATES *
 **********************/

// usage/help menu
void print_usage();

// output buffer functions
struct output_buffer* create_output_buffer(int size);
void output_buffer_append(struct output_buffer* output_buffer, char* buffer);
char* output_buffer_get(struct output_buffer* output_buffer, int i);
void free_output_buffer(struct output_buffer* output_buffer);

// output funtions
void print_motif_counts(int N_min, int N_max, int L_min, int L_max);
void print_all_motifs(int N_min, int N_max, int L_min, int L_max);
void print_distance_statistics(char* pdb_file, char* stems_file, char* suffix, double bin_size);
void print_C_eff_for_all_motifs(char* pdb_file, char* stems_file, char* suffix, int N_min, int N_max, int L_min, int L_max, int N_samples, int N_threads, double bin_size, double d_min, int self_avoiding, int sample_count, char* helix_file);
void print_C_eff_for_motif_file(char* pdb_file, char* stems_file, char* suffix, char* motif_file, int N_samples, int N_threads, double bin_size, double d_min, int self_avoiding, char* helix_file);
void print_C_eff_for_input_motif(char* pdb_file, char* stems_file, char* suffix, char* input_motif, int N_samples, int N_threads, double bin_size, double d_min, int self_avoiding, char* helix_file);
void print_loop_trace_for_input_motif(char* pdb_file, char* stems_file, char* suffix, char* input_motif, int N_samples,int N_threads, double bin_size, double d_min, int self_avoiding, char* helix_file, double p_circ_min);
/********
 * MAIN *
 ********/

int main(int argc, char *argv[]) {

  // default values
  double bin_size = 3.0; // should probably just force the bin size to be same as d_min
  double d_min = 3.0;
  double p_circ_min = 1e-4;
  int N_min = 2;
  int N_max = 5;
  int L_min = 0;
  int L_max = 4;
  int N_samples = 2000;
  int N_threads = 1;
  int sample_count = 0;
  int export_coords = 0;  
  int counts_only = 0;
  int dist_stats = 0;
  int self_avoiding = 0;
  char* input_motif = NULL;
  char* pdb_file = NULL;
  char* stems_file = NULL;
  char* suffix = NULL;
  char* motif_file = NULL;
  char* helix_file = NULL;
  char opt;
  
  while((opt = getopt(argc,argv,"n:r:N:l:L:b:c:d:m:s:t:x:P:S:M:H:ECDAh")) != -1) {
    switch(opt) {
      // First the run-modes
    case 'C':
      counts_only = 1;
      break;
    case 'D':
      dist_stats = 1;
      break;
    case 'A':
      self_avoiding = 1;
      break;
    case 'E':
      export_coords = 1;
      break;
    case 'n': // N_min, minimum number of loops to consider
      N_min = atoi(optarg);
      break;
    case 'N': // N_max, maximum number of loops to consider
      N_max = atoi(optarg);
      break;
    case 'l': // L_min the minimum number of loops to consider
      L_min = atoi(optarg);
      break;
    case 'L':  // L_max the maximum number of loops to consider
      L_max = atoi(optarg);
      break;
    case 's': // N_samples, the number of iterations/number of random trajectories to sample 
      N_samples = atoi(optarg);
      break;
    case 'r': // sample_count, the number of randomly sampled loops to sample if it exceeds this number
      sample_count = atoi(optarg);
      break;
    case 't': // N_threads, the number of threads used in OMP parallel programmaing over motifs
      N_threads = atoi(optarg);
      break;
    case 'b': // bin_size, the size of bins for the collision detector.
      bin_size = atof(optarg);
      break;
    case 'c': // minimum circularlization probability for exporting loop traces
      p_circ_min = atof(optarg);
      break;
    case 'd': // d_min, the minimum distance for collision detection
      d_min = atof(optarg);
      break;
      // The input files that can be provided.
    case 'm': // motif string
      input_motif = malloc((strlen(optarg)+1) * sizeof(char));
      strcpy(input_motif,optarg);
      break;
    case 'P': // PDB file, required for sampling transforms
      pdb_file = malloc((strlen(optarg)+1) * sizeof(char));
      strcpy(pdb_file,optarg);
      break;
    case 'S': // stems file, required for sampling transforms
      stems_file = malloc((strlen(optarg)+1) * sizeof(char));
      strcpy(stems_file,optarg);
      break;      
    case 'M': // motif file, required for sampling transforms
      motif_file = malloc((strlen(optarg)+1) * sizeof(char));
      strcpy(motif_file,optarg);
      break;
    case 'H': // helix file for creating idealized helix
      helix_file = malloc((strlen(optarg)+1) * sizeof(char));
      strcpy(helix_file,optarg);
      break;
    case 'x': // suffix to add to the output files
      suffix = malloc((strlen(optarg)+1) * sizeof(char));
      strcpy(suffix,optarg);
      break;
    case 'h':
      print_usage();
      break;
    default:
      printf("Unknown option: '%c'.\n",opt);
      print_usage();
      break;
    }
  }

  srand((unsigned int) time(NULL)); // set random seed by time
  //srand(7);
  int r = rand() % 7; // discard first value according to some stackoverflow post I saw
  
  if((pdb_file)&&(stems_file)) {
    if(helix_file && !self_avoiding) {
      printf("\n\nYou've provided a helix file, but this is for self-avoiding polymer mode only ( -a ).\n");
      print_usage();
    }
    if(motif_file) {      
      print_C_eff_for_motif_file(pdb_file,stems_file,suffix,motif_file,N_samples,N_threads,
				 bin_size,d_min,self_avoiding,helix_file);
    } else if(input_motif) {
      if(export_coords) {
        print_loop_trace_for_input_motif(pdb_file,stems_file,suffix,input_motif,N_samples,N_threads,bin_size,d_min,self_avoiding,helix_file,p_circ_min);
      } else {
        print_C_eff_for_input_motif(pdb_file,stems_file,suffix,input_motif,N_samples,N_threads,bin_size,d_min,self_avoiding,helix_file);
      }
    } else {
      if(dist_stats) {
	print_distance_statistics(pdb_file,stems_file,suffix,bin_size);
      } else {
	print_C_eff_for_all_motifs(pdb_file,stems_file,suffix,N_min,N_max,L_min,L_max,N_samples,N_threads,
				   bin_size,d_min,self_avoiding,sample_count,helix_file);
      }
    }
  } else if(motif_file && !(pdb_file && stems_file)) {
    printf("\n\nPlease enter a PDB file and stems file when inputing a motif file to calculate C_eff values.\n");
    print_usage();
  } else if(input_motif && !(pdb_file && stems_file)) {
    printf("\n\nPlease enter a PDB file and stems file when inputing a motif to calculate C_eff values.\n");
    print_usage();
  } else if(helix_file && !(pdb_file && stems_file)) {
    printf("\n\nYou've provided a helix file, but his also requires a main PDB file and stems file.\n");
    print_usage();
  } else if(helix_file && !self_avoiding) {
    printf("\n\nYou've provided a helix file, but this is for self-avoiding polymer mode only ( -a ).\n");
    print_usage();    
  } else if(counts_only) {
    print_motif_counts(N_min,N_max,L_min,L_max);
  } else {
    print_all_motifs(N_min,N_max,L_min,L_max);
  }
  if(pdb_file) free(pdb_file);
  if(stems_file) free(stems_file);
  if(motif_file) free(motif_file);
  if(helix_file) free(helix_file);
  if(input_motif) free(input_motif);
  if(suffix) free(suffix);
  return 0;
}

/********************
 * Usage/Help Menu  *
 ********************/

void print_usage() {
  printf("\nUsage: \n\n");
  printf("./paperloops -P <RNA PDB_file> -S <stems file> -n <N_min> -N <N_max> -l <L_min> -L <L_max> [options]\n\n");
  printf("Options:\n");
  printf("\t-C : Enumerate the number of possible loop motifs and print to a file.\n\n");
  printf("\t-D : Print distance statistics only for the provided PDB file\n\n");  
  printf("\t-A : Run in self-avoiding model (\"real\" polymers, not ideal).\n\n");
  printf("\t-E : Export coordinates of circularizing RNA loops to file [loop_traces_<suffix>.txt]\n\n");
  printf("\t-n <int> : N_min, the minimum number of loops to consider. [default: 2]\n\n");
  printf("\t-N <int> : N_max, the maximum number of loops to consider. [default: 5]\n\n");
  printf("\t-l <int> : L_min, the minimum total loop length to consider. [default: 0]\n\n");
  printf("\t-L <int> : L_max, the maximum total loop length to consider. [default: 4]\n\n");
  printf("\t-s <int> : N_samples, the number of iterations/random trajectories to consider [default: 2000].\n\n");
  printf("\t-r <int> : sample_count, the number of random samples to select among possible loops for n,L combinations. [default: 0, no sampling]\n\n");
  printf("\t-t <int> : N_threads, the number of threads to use for parallel implementation [default 1]\n\n");
  printf("\t-b <float> : bin_size, the bin size for the detector patches.\n\n");
  printf("\t-c <float> : p_circ_min, the minimum circularization probability to export to file (used with -E) [default: 1e-4].\n\n");
  printf("\t-d <float> : d_min, the minimum distance to be considered a collision [default: 3.0 Angstroms].\n\n");
  printf("\t-m <string> : input_motif, a single input motif to be used for a solo paperloops run (e.g. 3-2-1).\n\n");
  printf("\t-P <PDB file> : the path to a PDB file to read in. [default: None]\n\n");
  printf("\t-S <Stems file> : the path to a stems file generated by rna_motif for PDB file. [default: None]\n\n");
  printf("\t-M <Motif file> : the path to a motif file to compute C_eff for instead of enumerating all of them.\n\n");
  printf("\t-H <Helix file> : the path to a PDB file with an ideal helix: 1rna.pdb.\n\n");
  printf("\t-x <output suffix> : a suffix added to the output files. [default: No suffix added]\n\n");
  printf("\t-h : print this help menu\n\n");
  exit(1);
}

/*******************
 * OUTPUT BUFFER   *
 *******************/

struct output_buffer* create_output_buffer(int size) {
  struct output_buffer* output_buffer = malloc(sizeof(struct output_buffer));
  output_buffer->len = size;
  output_buffer->count = 0;
  output_buffer->data = malloc(size * sizeof(char*));
  return output_buffer;
}

void free_output_buffer(struct output_buffer* output_buffer) {
  assert(output_buffer);
  assert(output_buffer->data);
  // free the two places where malloc was run:
  int count = output_buffer->count;
  for(int i=0; i<count; i++) {
    free(output_buffer->data[i]);
  }
  free(output_buffer->data);
  free(output_buffer);
}

void output_buffer_append(struct output_buffer* output_buffer, char* buffer) {
  assert(output_buffer);
  int i = output_buffer->count;
  output_buffer->data[i] = malloc((strlen(buffer)+1)*sizeof(char));
  strcpy(output_buffer->data[i], buffer);
  output_buffer->count++;
  //for(int j=0; j<=i; j++) {
  //  printf("%d %s\n",j,output_buffer->data[j]);
  //}
}

void output_buffer_set(struct output_buffer* output_buffer, int i, char* buffer) {
  assert(output_buffer);
  output_buffer->data[i] = malloc((strlen(buffer)+1)*sizeof(char));
  strcpy(output_buffer->data[i], buffer);
  output_buffer->count++;
}

char* output_buffer_get(struct output_buffer* output_buffer, int i) {
  assert(output_buffer);
  int count = output_buffer->count;
  int len = output_buffer->len;
  assert((0 <= i)&&(i < count)&&(i < len));
  return output_buffer->data[i];
}

/********************
 * OUTPUT FUNCTIONS *
 ********************/

void print_motif_counts(int N_min, int N_max, int L_min, int L_max) {
  // print an output file of all theoretical motif counts based on equation
  printf("Enumerating possible motifs for N_min=%d to N_max=%d and L_min=%d to L_max=%d.\n",N_min,N_max,L_min,L_max);
  char* out_file = "motif_counts.txt";
  FILE* OUT = fopen(out_file,"w");
  printf("printing out motif counts to: %s.\n",out_file);

  int** N = necklace_counts_fixed_length(N_max,L_max);
  fprintf(OUT,"n\tL\tN(n,L)\n");
  for(int n=N_min; n<=N_max; n++) {
    for(int L=L_min; L<=L_max; L++) {
      fprintf(OUT,"%d\t%d\t%d\n",n,L,N[n][L]);
    }
  }
  free_necklace_counts(N,N_max);
  fclose(OUT);
}

void print_all_motifs(int N_min, int N_max, int L_min, int L_max) {
  // print an output file of all motifs for input ranges
  // theoretical counts based on equation
  printf("Enumerating possible motifs for N_min=%d to N_max=%d and L_min=%d to L_max=%d.\n",N_min,N_max,L_min,L_max);
  
  int** N = necklace_counts_fixed_length(N_max,L_max);
  for(int n=N_min; n<=N_max; n++) {
    for(int L=L_min; L<=L_max; L++) {
      printf("# Npred(n=%d,L=%d)=%d\n",n,L,N[n][L]);
    }
  }
  // actual motifs built with recurssive efficient algorithm
  printf("Building motif file list.\n");
  struct motif_group*** M = build_motif_lib(N_min,N_max,L_max,N);
  
  char* out_file = "loop_motif_list.txt";
  FILE* OUT = fopen(out_file,"w");
  printf("printing out motifs to: %s.\n",out_file);

  for(int n=N_min; n<=N_max; n++) {
    for(int L=L_min; L<=L_max; L++) {
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
  free_motif_lib(M,N_max,L_max);
  free_necklace_counts(N,N_max);
  fclose(OUT);
}

void print_distance_statistics(char* pdb_file, char* stems_file, char* suffix, double bin_size) {
  printf("reading in PDB file: %s\n",pdb_file);
  struct rna_data* RNA = read_pdb_file(pdb_file);
  
  // get the closest approach:
  double d_min_obs = get_rna_closest_approach(RNA,bin_size);
  printf("closest approach: %.2f\n",d_min_obs);
}

/*********************
 * CORE CALCULATIONS *
 *********************/

void print_C_eff_for_all_motifs(char* pdb_file, char* stems_file, char* suffix, int N_min, int N_max, int L_min, int L_max, int N_samples, int N_threads, double bin_size, double d_min, int self_avoiding, int sample_count, char* helix_file) {

  // theoretical motif counts
  printf("Enumerating possible motifs for N_min=%d to N_max=%d and L_min=%d to L_max=%d.\n",N_min,N_max,L_min,L_max);
  int** N = necklace_counts_fixed_length(N_max,L_max);

  // actual motifs built with efficient algorithm
  printf("Building motif file list.\n");
  struct motif_group*** M = build_motif_lib(N_min,N_max,L_max,N);
  
  // building libraries
  struct library* rna_lib = build_library(pdb_file, stems_file);
  // building helix libraries
  struct library* helix_lib = build_helix_library(helix_file);
  
  // define output files, account for possible suffix
  char out_all_file[50] = "C_eff_all_results";
  char out_file[50] = "C_eff_results";
  char sap_report[50] = "collision_report"; // defined even for non-SAP mode
  if(suffix) {
    strcat(out_all_file,"_");
    strcat(out_file,"_");
    strcat(sap_report,"_");
    strcat(out_all_file,suffix);
    strcat(out_file,suffix);
    strcat(sap_report,suffix);
  }
  strcat(out_all_file,".txt");
  strcat(out_file,".txt");
  strcat(sap_report,".txt");

  FILE* OUT_ALL = fopen(out_all_file,"w");
  FILE* OUT = fopen(out_file,"w");
  if(self_avoiding) {
    // clear the sap_report file if about to append to it 
    FILE* SAP_REPORT = fopen(sap_report,"w");
    fprintf(SAP_REPORT,"motif\tsteptypes\tprob\tpair_collide\tpair_tot\tfwd_collide\trev_collide\tN_samples\n");
    fclose(SAP_REPORT);
  }
  printf("Writing C_eff results to: %s and %s\n",out_file,out_all_file);

  // inclusive loops over the provided boundaries, motifs already built
  for(int n=N_min; n<=N_max; n++) {
    for(int L=L_min; L<=L_max; L++) {
      // update standard output log
      int count = get_motif_group_count(M[n][L]);
      printf("# N_pred(n=%d,L=%d)=%d, N_obs=%d\n",n,L,N[n][L],count);
      
      // sample the motifs if needed
      if(sample_count && count > sample_count) {
        // Knuth shuffle - partial Fisher-Yates shuffle
        for(int k=0; k<sample_count; k++) {
          // random sample from k to count-1
          int range = count - k;
          int random_index = k + rand() % range;
          struct motif* m1 = get_motif(M[n][L],k);
          struct motif* m2 = get_motif(M[n][L],random_index);
          set_motif(M[n][L],k,m2);
          set_motif(M[n][L],random_index,m1);            
        }
        // update the count to just get sample_count examples
        count = sample_count;
      }
      
      for(int i=0; i<count; i++) {
	// here we have a motif, part of a necklace, first one should be lexicographically lowest.
        // attached to it are all distinct rotations in a linked list rot1->rot2->rot3->...
	struct motif* m = get_motif(M[n][L],i);
	// loop through the motif rotations
	while(m) {
	  struct steptypes* steps = motif_to_steptypes(m);
	  int steps_length = steptypes_length(steps);
	  // prepare the output buffer
          //printf("allocating %d sized buffer\n",steps_length);
	  struct output_buffer* output_buffer = create_output_buffer(steps_length);
	  // initialize for averages
	  double num = 0.0;
	  double den = 0.0;
          double full_dist_sq_mean = 0.0;
          double full_theta_sq_mean = 0.0;

	  // average over the rotations of the steptypes:
          #pragma omp parallel for num_threads(N_threads) default(none) \
          shared(i, n, L, steps, steps_length, m, rna_lib, num, den, full_dist_sq_mean, full_theta_sq_mean, N_samples, bin_size, d_min, helix_lib, sap_report, self_avoiding, output_buffer) \
          schedule(dynamic, 1)
	  for(int j=0; j<steps_length; j++) {

            // thread-safe copy of steptypes
            struct steptypes* thread_steps = clone_steptypes(steps); 
            for(int r=0; r<j; r++) {
              thread_steps = rotate_steptypes(thread_steps);
            }

            //printf("in (%d): ",j);
            //print_steptypes(steps);
            //printf("out(%d): ",j);
            //print_steptypes(thread_steps);
            
            // allocate tracking struct on the stack
	    double C_eff_stats[4];
            
	    if(self_avoiding) {
	      get_C_eff_both_dir_sap(m,thread_steps,rna_lib,N_samples,bin_size,d_min,helix_lib,sap_report,C_eff_stats);
	    } else {
	      get_C_eff_both_dir(thread_steps,rna_lib,N_samples,C_eff_stats);
	    }
            // unpack the stats
	    double C_eff_mean = C_eff_stats[0];
	    double C_eff_error = C_eff_stats[1];
            double dist_sq_mean = C_eff_stats[2];
            double theta_sq_mean = C_eff_stats[3];
            
	    // write stats to the output buffer
	    char buffer[500]; // stack buffer to print the stats for this record
            int offset = sprintf(buffer,"%d\t%d\t%d\t",i,n,L);
	    offset += sprint_motif(buffer+offset,m,'\t');
	    offset += sprint_steptypes(buffer+offset,thread_steps,'\t');
	    offset += sprintf(buffer+offset,"%.8f\t%.8f\t%.8f\t%.8f\n",C_eff_mean,C_eff_error,dist_sq_mean,theta_sq_mean);

            #pragma omp critical
            {
              output_buffer_set(output_buffer,j,buffer);
              // compute relative error 
              double C_eff_relerr = C_eff_error / C_eff_mean;	    
              // logarithmic, inverse-relerr weighted average
              num += log(C_eff_mean)/(C_eff_relerr*C_eff_relerr);
              den += 1.0/(C_eff_relerr*C_eff_relerr);
              // averages over step-type rotations
              full_dist_sq_mean += dist_sq_mean;
              full_theta_sq_mean += theta_sq_mean;
            }
            // memory management:
            free_steptypes(thread_steps);
	  }
	  // weighted average over steptypes rotations:
	  double full_C_eff_mean = exp(num/den);
	  double full_C_eff_err = sqrt(1.0/den)*full_C_eff_mean;
          full_dist_sq_mean /= steps_length;
          full_theta_sq_mean /= steps_length;
          
          // print all motif steptype rotations, result stored in output_buffer
          for(int j=0; j<steps_length; j++) {
            fprintf(OUT_ALL, "%s", output_buffer_get(output_buffer,j));
          }

          // print average of all rotations
          fprintf(OUT, "%d\t%d\t%d\t",i,n,L);
          fprint_motif(OUT, m, '\t');
          fprintf(OUT,"%.8f\t%.8f\t%.8f\t%.8f\n",full_C_eff_mean,full_C_eff_err,full_dist_sq_mean,full_theta_sq_mean);    
          // rotate the motif to get an equivalent one in different rotation
          m = get_next_motif(m);
          // memory management
          free_steptypes(steps);
          free_output_buffer(output_buffer);
        }
      }
    }
  }
  free_necklace_counts(N,N_max);
  free_motif_lib(M,N_max,L_max);
  free_library(rna_lib);
  free_library(helix_lib);
  fclose(OUT);
  fclose(OUT_ALL);
}

void print_C_eff_for_input_motif(char* pdb_file, char* stems_file, char* suffix, char* input_motif, int N_samples, int N_threads, double bin_size, double d_min, int self_avoiding, char* helix_file) {

  // building libraries
  struct library* rna_lib = build_library(pdb_file, stems_file);
  // building helix libraries
  struct library* helix_lib = build_helix_library(helix_file);
  
  // define output files, account for possible suffix
  char out_all_file[50] = "C_eff_all_results.txt";
  char out_file[50] = "C_eff_results.txt";
  char sap_report[50] = "collision_report.txt"; // defined even for non SAP mode
  if(suffix) {
    strcpy(out_all_file,"C_eff_all_results_");
    strcpy(out_file,"C_eff_results_");
    strcpy(sap_report,"collision_report_");
    strcat(out_all_file,suffix);
    strcat(out_file,suffix);
    strcat(sap_report,suffix);
    strcat(out_all_file,".txt");
    strcat(out_file,".txt");
    strcat(sap_report,".txt");
  }
  FILE* OUT_ALL = fopen(out_all_file,"w");
  FILE* OUT = fopen(out_file,"w");
  if(self_avoiding) {
    // clear the sap_report file if about to append to it 
    FILE* SAP_REPORT = fopen(sap_report,"w");
    fprintf(SAP_REPORT,"motif\tsteptypes\tprob\tpair_collide\tpair_tot\tfwd_collide\trev_collide\tN_samples\n");
    fclose(SAP_REPORT);
  }
  printf("Writing C_eff results to: %s and %s\n",out_file,out_all_file);

  printf("Reading input motif: %s\n",input_motif);
  struct motif* m = read_input_motif(input_motif);  
  //print_motif(m,'\n');
    
  // loop through the motif rotations
  while(m) {
    struct steptypes* steps = motif_to_steptypes(m);	  	  
    int steps_length = steptypes_length(steps);
    // prepare the output buffer
    struct output_buffer* output_buffer = create_output_buffer(steps_length);
    //initialize for averages
    double num = 0.0;
    double den = 0.0;
    double full_dist_sq_mean = 0.0;
    double full_theta_sq_mean = 0.0;
    
    // average over the rotations of the steptypes:
    for(int j=0; j<steps_length; j++) {
      //printf("C_eff for motif: ");
      //print_motif(m,'\n');
      double C_eff_stats[4];
      
      if(self_avoiding) {
	  get_C_eff_both_dir_sap(m,steps,rna_lib,N_samples,bin_size,d_min,helix_lib,sap_report,C_eff_stats);
      } else {
        get_C_eff_both_dir(steps,rna_lib,N_samples,C_eff_stats);
      }
      double C_eff_mean = C_eff_stats[0];
      double C_eff_error = C_eff_stats[1];
      double dist_sq_mean = C_eff_stats[2];
      double theta_sq_mean = C_eff_stats[3];
      
      // write stats to the output buffer
      char buffer[300]; // allocated buffer to print the stats for this record
      int offset = sprint_motif(buffer,m,'\t');
      offset += sprint_steptypes(buffer+offset,steps,'\t');
      offset += sprintf(buffer+offset,"%.8f\t%.8f\t%.8f\t%.8f\n",C_eff_mean,C_eff_error,dist_sq_mean,theta_sq_mean);
      output_buffer_set(output_buffer,j,buffer);
      
      // averages over step-type rotations
      full_dist_sq_mean += dist_sq_mean;
      full_theta_sq_mean += theta_sq_mean;
      // compute relative error 
      double C_eff_relerr = C_eff_error / C_eff_mean;
      // logarithmic, inverse-relerr weighted average
      num += log(C_eff_mean)/(C_eff_relerr*C_eff_relerr);
      den += 1.0/(C_eff_relerr*C_eff_relerr);
      // rotate the steptypes and add to running average
      steps = rotate_steptypes(steps);
    }
    // weighted average over steptypes rotations:
    double full_C_eff_mean = exp(num/den);
    double full_C_eff_err = sqrt(1.0/den)*full_C_eff_mean;
    full_dist_sq_mean /= steps_length;
    full_theta_sq_mean /= steps_length;
    
    for(int j=0; j<steps_length; j++) {
      fprintf(OUT_ALL, "%s", output_buffer_get(output_buffer,j));
    }
    free_output_buffer(output_buffer);
    fprint_motif(OUT,m,'\t');
    fprintf(OUT,"%.8f\t%.8f\t%.8f\t%.8f\n",full_C_eff_mean,full_C_eff_err,full_dist_sq_mean,full_theta_sq_mean);

    m = get_next_motif(m); // rotate the motif to get an equivalent one in different rotation    
    free_steptypes(steps);      
  }
  free_library(rna_lib);
  free_library(helix_lib);
  fclose(OUT);
  fclose(OUT_ALL);
}

void print_C_eff_for_motif_file(char* pdb_file, char* stems_file, char* suffix, char* motif_file, int N_samples, int N_threads, double bin_size, double d_min, int self_avoiding, char* helix_file) {

  // building libraries
  struct library* rna_lib = build_library(pdb_file, stems_file);
  // building helix libraries
  struct library* helix_lib = build_helix_library(helix_file);
  
  printf("Reading motif file: %s\n",motif_file);
  struct motif_group* motif_group = read_motif_file(motif_file);

  // define output files, account for possible suffix
  char out_all_file[50] = "C_eff_all_results.txt";
  char out_file[50] = "C_eff_results.txt";
  char sap_report[50] = "collision_report.txt"; // defined even for non SAP mode
  if(suffix) {
    strcpy(out_all_file,"C_eff_all_results_");
    strcpy(out_file,"C_eff_results_");
    strcpy(sap_report,"collision_report_");
    strcat(out_all_file,suffix);
    strcat(out_file,suffix);
    strcat(sap_report,suffix);
    strcat(out_all_file,".txt");
    strcat(out_file,".txt");
    strcat(sap_report,".txt");
  }
  FILE* OUT_ALL = fopen(out_all_file,"w");
  FILE* OUT = fopen(out_file,"w");
  if(self_avoiding) {
    // clear the sap_report file if about to append to it 
    FILE* SAP_REPORT = fopen(sap_report,"w");
    fprintf(SAP_REPORT,"motif\tsteptypes\tprob\tpair_collide\tpair_tot\tfwd_collide\trev_collide\tN_samples\n");
    fclose(SAP_REPORT);
  }
  printf("Writing C_eff results to: %s and %s\n",out_file,out_all_file);

  int count = get_motif_group_count(motif_group);
  printf("count = %d motifs\n",count);

  for(int i=0; i<count; i++) {
    // extract the i-th motif
    struct motif* m = get_motif(motif_group,i);
    // loop through the motif rotations
    while(m) {
      struct steptypes* steps = motif_to_steptypes(m);	  	  
      int steps_length = steptypes_length(steps);
      // prepare the output buffer
      struct output_buffer* output_buffer = create_output_buffer(steps_length);
      //initialize for averages
      double num = 0.0;
      double den = 0.0;
      double full_dist_sq_mean = 0.0;
      double full_theta_sq_mean = 0.0;
      
      // average over the rotations of the steptypes:
      for(int j=0; j<steps_length; j++) {
	double C_eff_stats[4];

        if(self_avoiding) {
	  get_C_eff_both_dir_sap(m,steps,rna_lib,N_samples,bin_size,d_min,helix_lib,sap_report,C_eff_stats);
	} else {
	  get_C_eff_both_dir(steps,rna_lib,N_samples,C_eff_stats);
	}
	double C_eff_mean = C_eff_stats[0];
	double C_eff_error = C_eff_stats[1];
        double dist_sq_mean = C_eff_stats[2];
        double theta_sq_mean = C_eff_stats[3];

	// write stats to the output buffer
	char buffer[300]; // allocated buffer to print the stats for this record
	int offset = sprint_motif(buffer,m,'\t');
	offset += sprint_steptypes(buffer+offset,steps,'\t');
        offset += sprintf(buffer+offset,"%.8f\t%.8f\t%.8f\t%.8f\n",C_eff_mean,C_eff_error,dist_sq_mean,theta_sq_mean);
	output_buffer_set(output_buffer,j,buffer);
	
        // averages over step-type rotations
        full_dist_sq_mean += dist_sq_mean;
        full_theta_sq_mean += theta_sq_mean;
	// compute relative error 
	double C_eff_relerr = C_eff_error / C_eff_mean;
	// logarithmic, inverse-relerr weighted average
	num += log(C_eff_mean)/(C_eff_relerr*C_eff_relerr);
	den += 1.0/(C_eff_relerr*C_eff_relerr);
	// rotate the steptypes and add to running average
	steps = rotate_steptypes(steps);
	// memory management:
      }
      // weighted average over steptypes rotations:
      double full_C_eff_mean = exp(num/den);
      double full_C_eff_err = sqrt(1.0/den)*full_C_eff_mean;
      full_dist_sq_mean /= steps_length;
      full_theta_sq_mean /= steps_length;

      for(int j=0; j<steps_length; j++) {
        fprintf(OUT_ALL, "%s", output_buffer_get(output_buffer,j));
      }
      free_output_buffer(output_buffer);
      fprint_motif(OUT,m,'\t');
      fprintf(OUT,"%.8f\t%.8f\t%.8f\t%.8f\n",full_C_eff_mean,full_C_eff_err,full_dist_sq_mean,full_theta_sq_mean);

        m = get_next_motif(m); // rotate the motif to get an equivalent one in different rotation
      free_steptypes(steps);      
    }
  }
  free_motif_group(motif_group);
  free_library(rna_lib);
  free_library(helix_lib);
  fclose(OUT);
  fclose(OUT_ALL);
}

void print_loop_trace_for_input_motif(char* pdb_file, char* stems_file, char* suffix, char* input_motif, int N_samples, int N_threads, double bin_size, double d_min, int self_avoiding, char* helix_file, double p_circ_min) {
  // building libraries
  struct library* rna_lib = build_library(pdb_file, stems_file);
  // building helix libraries
  struct library* helix_lib = build_helix_library(helix_file);
  
  // define output files, account for possible suffix
  char out_file[50] = "loop_traces.txt";
  char sap_report[50] = "collision_report.txt"; // defined even for non SAP mode
  if(suffix) {
    strcpy(out_file,"loop_traces_");
    strcpy(sap_report,"collision_report_");
    strcat(out_file,suffix);
    strcat(sap_report,suffix);
    strcat(out_file,".txt");
    strcat(sap_report,".txt");
  }
  FILE* OUT = fopen(out_file,"w");
  if(self_avoiding) {
    // clear the sap_report file if about to append to it 
    FILE* SAP_REPORT = fopen(sap_report,"w");
    fprintf(SAP_REPORT,"motif\tsteptypes\tprob\tpair_collide\tpair_tot\tfwd_collide\trev_collide\tN_samples\n");
    fclose(SAP_REPORT);
  }
  printf("Writing loop traces to: %s\n",out_file);

  printf("Reading input motif: %s\n",input_motif);
  struct motif* m = read_input_motif(input_motif);  
    
  // loop through the motif rotations
  int idx = 0;
  while(m) {
    struct steptypes* steps = motif_to_steptypes(m);	  	  
    int steps_length = steptypes_length(steps);
    
    // average over the rotations of the steptypes:
    for(int j=0; j<steps_length; j++) {
      if(self_avoiding) {
        print_loop_traces_sap(OUT,&idx,m,steps,rna_lib,N_samples,bin_size,d_min,helix_lib,sap_report,p_circ_min);
      } else {
        print_loop_traces(OUT,&idx,m,steps,rna_lib,N_samples,helix_lib,p_circ_min);
      }
      // rotate the steptypes
      steps = rotate_steptypes(steps);
    }
    m = get_next_motif(m); // rotate the motif to get an equivalent one in different rotation    
    free_steptypes(steps);      
  }
  free_library(rna_lib);
  free_library(helix_lib);
  fclose(OUT);
}
