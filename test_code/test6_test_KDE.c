#include "rna_io.h"
#include "rna_data.h"
#include "simulation.h"
#include "transform.h"
#include "motifs.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

//struct library* read_pts_file(char* pts_file);

int main(int argc, char *argv[]) {
  if(argc == 4) {
    ;
  } else {
    printf("Usage: %s <pts fwd file> <pts rev file> <suffix>\n",argv[0]);
    exit(1);
  }
  char* pts_fwd_file = argv[1];
  char* pts_rev_file = argv[2];
  char* suffix = argv[3];
  
  srand(time(NULL)); // set random seed by time
  //srand(777);
  int r = rand() % 7; // discard first value based on stackoverflow
  
  struct library* pts_fwd = read_pts_file(pts_fwd_file);
  struct library* pts_rev = read_pts_file(pts_rev_file);
  struct array* h_fwd = get_kde_bandwidth(pts_fwd);
  
  FILE* PROBS = fopen("test_probs.txt","w");
  for(int j=0; j<6; j++) {
    double h = array_get(h_fwd,j);
    fprintf(PROBS,"bandwidth=%.3f\n",h);
  }
  int N = library_size(pts_rev);
  for(int i=0; i<N; i++) {
    struct transform* x = get_transform(pts_rev,i);
    double p = get_kde_probability(x, pts_fwd, h_fwd);
    fprintf(PROBS,"p=%.6e\n",p);
  }
  fclose(PROBS);
  
  return 0;
}

/*
struct library* read_pts_file(char* pts_file) {
  int MAX_LINE_LENGTH = 1000;
  FILE* F = fopen(pts_file,"r");
  if(F == NULL) {
    printf("Error opening the file:\"%s\".\n Check permissions/file name.\n",pts_file);
    exit(1);
  }

  struct library* pts = create_library();
  char line[MAX_LINE_LENGTH];
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    struct split_text* row = split_line(line,"\t");
    if(get_split_text_cols(row) == 6) {
      struct transform* T = create_transform();
      double* t = create_vector();
      set_vector_elements(t,atof(get_split_text_term(row,0)),atof(get_split_text_term(row,1)),atof(get_split_text_term(row,2)));
      set_translation_vector(T,t);
      double* EV = create_vector();
      set_vector_elements(EV,atof(get_split_text_term(row,3)),atof(get_split_text_term(row,4)),atof(get_split_text_term(row,5)));
      set_vector(get_euler_vector(T),EV);
      append_transform(pts,T); 
    }
  }
  return pts;
}
*/
