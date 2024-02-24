/*
 * This library contains the definitions of structures and functions reading PDB files
 * for RNA structures and storing into the RNA data structure defined in RNA_molecule.c
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "motifs.h"
#include "rna_io.h"
#include "rna_data.h"

#define MAX_LINE_LENGTH 1000

/*******************
 * data structures *
 *******************/

struct split_text {
  char** term;
  int cols;
};

/* FUNCTIONS */

double* extract_coord(char* line) {
  // extract the substrings for x,y,z
  char x_str[9];
  sprintf(x_str,"%.8s",line+30);
  char y_str[9];
  sprintf(y_str,"%.8s",line+38);
  char z_str[9];
  sprintf(z_str,"%.8s",line+46);
  // create an array of doubles to return the coord
  double *r = malloc(3*sizeof(double));
  r[0] = atof(x_str);
  r[1] = atof(y_str);
  r[2] = atof(z_str);
  return r;
}

int extract_atom_num(char* line) {
  // should be pos 7-11 in 1-based, 6-10 0-based
  char atom_num_str[6];
  sprintf(atom_num_str,"%.5s",line+6);
  return atoi(atom_num_str);
}

char* extract_atom_name(char* line) {
  // should be pos 13-16 in 1-based
  char* atom_name = malloc(5*sizeof(char));
  int count = 0;
  for(int i=0;i<4;i++) {
    if(line[12+i] != ' ') {
      atom_name[count++] = line[12+i];
    }
  }
  atom_name[count] = '\0';
  //sprintf(atom_name,"%.4s",line+12);
  return atom_name;
}

int extract_nuc_num(char* line) {
  // should be pos 23-26 in 1-based
  char nuc_num_str[5];
  sprintf(nuc_num_str,"%.4s",line+22);
  return atoi(nuc_num_str);
}

char* extract_nuc_name(char* line) {
  // should be pos 18-20 in 1-based
  char* nuc_name = malloc(4*sizeof(char));
  int count = 0;
  for(int i=0;i<3;i++) {
    if(line[17+i] != ' ') {
      nuc_name[count++] = line[17+i];
    }
  }
  nuc_name[count] = '\0';
  // NEED to map PSU to P and other HETATM nucleotides
  return nuc_name;
}

char extract_chain_id(char* line) {
  char chain_id = line[21];
  return chain_id;
}

int read_model_count(FILE* F) {
  char line[MAX_LINE_LENGTH];
  
  int model_count = 1; // one model by default 
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    if(strncmp(line,"MODEL ",6) == 0) {
      // MODEL        1                                                       
      char model_num_str[5];
      sprintf(model_num_str,"%.4s",line+11); // model serial number at positions 11-14
      int model_num = atoi(model_num_str);
      // here I am assuming no model numbers are "skipped". 
      if(model_num > model_count) {
	model_count = model_num;
      }
    }
  }
  printf("Found %d models.\n",model_count);
  return model_count;
}

int* read_model_locs(FILE* F, struct rna_data* RNA) {
  char line[MAX_LINE_LENGTH];
  int FOUND_MODEL = 0;
  int FOUND_ATOM = 0;

  // get the number of models 
  int model_count = read_model_count(F);
  set_model_count(RNA,model_count);
  
  // reset file and allocate memory
  fseek(F,0,0);
  allocate_model_locs(RNA);

  // read the model locations
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    if(strncmp(line,"ATOM  ",6) == 0) {
      if(!FOUND_MODEL && !FOUND_ATOM) {
	FOUND_ATOM = 1;
	int model_num = 0;
	int offset = ftell(F) - strlen(line);
	set_model_loc(RNA,model_num,offset);
      }
    } else if(strncmp(line,"MODEL ",6) == 0) {
      // MODEL        1                                                                  
      FOUND_MODEL = 1;
      char model_num_str[5];
      // model serial number at positions 11-14
      sprintf(model_num_str,"%.4s",line+11);
      int model_num = atoi(model_num_str);
      int offset = ftell(F);
      set_model_loc(RNA,model_num-1,offset); // convert model_num to 0-based for loading
    }
  }
}

int read_chain_count(FILE* F, struct rna_data* RNA, int model_num) {
  char line[MAX_LINE_LENGTH];
  fseek(F,get_model_loc(RNA,model_num),0);  // jump to the first ATOM line of model
  char prev_chain_id = '?';
  int COUNT1 = 0;
  int COUNT2 = 0; // second method of counting: double checking
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    if(strncmp(line,"ATOM  ",6) == 0) {
      char chain_id = extract_chain_id(line);
      if(prev_chain_id == '?') {
	prev_chain_id = chain_id;
	COUNT1++; // add to count for first chain
      }
      if(prev_chain_id != chain_id) {
	COUNT1++; // add for all other chains
	prev_chain_id = chain_id;
      }
    } else if(strncmp(line,"TER   ",6)) {
      COUNT2++;
    } else if(strncmp(line,"ENDMDL",6) == 0) {
      break;
    }
  }
  if(COUNT1 != COUNT2) {
    printf("chain counts don't match: %d %d\n",COUNT1,COUNT2);
    exit(1);
  }
  return COUNT1;
}

struct rna_data* read_pdb_file(char* input_file_name) {
  char line[MAX_LINE_LENGTH];
  
  FILE* F = fopen(input_file_name,"r");
  if(F == NULL) {
    printf("Error opening the file:\"%s\".\n Check permissions/file name.\n",input_file_name);
    exit(1);
  }

  struct rna_data* RNA = create_rna_data(input_file_name);
  // get the number and locations of each model, allocate memory
  read_model_locs(F, RNA);
  allocate_rna_models(RNA);
  
  for(int model_num=0; model_num < get_model_count(RNA); model_num++) {
    //int chain_count = read_model_chain_count(FILE,RNA,model_num); // get the number of chains for this model
    struct rna_model* RNA_model = create_model(); // create RNA_model
    // loop through the model locations and read in data: 
    fseek(F,get_model_loc(RNA,model_num),0);  // jump to the first ATOM line of model
    while(fgets(line,MAX_LINE_LENGTH,F)) {
      if((strncmp(line,"ATOM  ",6) == 0)||(strncmp(line,"HETATM",6) == 0)) {
	char chain_id = extract_chain_id(line);
	if(new_chain_id(RNA_model,chain_id)) {
	  printf("Found new chain: %c\n",chain_id);
	  struct chain* chain = create_chain(chain_id);
	  add_chain(RNA_model,chain);
	  //printf("added new chain: %c\n",chain_id);
	}
	int nuc_num = extract_nuc_num(line);
	char* nuc_name = extract_nuc_name(line);
	if(new_nuc_num(RNA_model,nuc_num)) {
	  struct nuc* n = create_nuc(nuc_num,nuc_name);
	  add_nuc(RNA_model,n);	  
	  //printf("Adding: %c %d %s\n",chain_id,nuc_num,nuc_name);
	  //set_chain_sequence(get_last_chain(RNA_model));
	  //printf("Current: %s\n",get_chain_sequence(get_last_chain(RNA_model)));
	}
	int atom_num = extract_atom_num(line);
	char* atom_name = extract_atom_name(line);
	double* r = extract_coord(line);
	struct atom* a = create_atom(atom_num,atom_name,r);
	add_atom(RNA_model,a);
	//printf("added atom: %s\n",atom_name);
      } else if(strncmp(line,"TER   ",6) == 0) {
	// expecting end of chain
	;
      } else if(strncmp(line,"ENDMDL",6) == 0) {
	break; // end of model
      }
    }
    // check number of chains.
    set_rna_model(RNA,model_num,RNA_model);    
  }
  // make an array copy of the linked list for faster access of base pairs
  set_rna_nuc_array(RNA);
  fclose(F);
  printf("PDB file %s read completed.\n",input_file_name);
  return RNA;
}

/*****************************
 * base pair reading methods *
 *****************************/

struct split_text* create_split_text() {
  struct split_text* row = malloc(sizeof(struct split_text));
  row->term = NULL;
  row->cols = 0;
  return row;
}

void free_split_text(struct split_text* row) {
  free(row->term);
  free(row);
}

int get_split_text_cols(struct split_text* row) {
  assert(row);
  return row->cols;
}

char* get_split_text_term(struct split_text* row, int i) {
  assert(row);
  return row->term[i];
}

struct split_text* split_line(char *line, char *delim) {
  struct split_text* row = create_split_text();
  int n = 0;
  
  char *token = strtok(line, delim);
  while( token != NULL ) {
    row->term = realloc(row->term, (n+1)*sizeof(char*));
    row->term[n] = token; 
    token = strtok(NULL, delim);
    n++;
  }
  row->cols = n;
  return row;
}

int count_basepairs(FILE* F) {
  char line[MAX_LINE_LENGTH];
  int size = 0;
  
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    struct split_text* row = split_line(line," ");
    if(row->cols >= 4) {
      if((strcmp(row->term[2],"W")==0)&&(strcmp(row->term[3],"W")==0)&&(strcmp(row->term[4],"C")==0)) {
	size++;
      }
    }
    free_split_text(row);
  }
  return size;
}

int count_stem_basepairs(FILE* F) {
  char line[MAX_LINE_LENGTH];
  int size = 0;
  
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    struct split_text* row = split_line(line," ");
    if(row->cols == 2) {
      //printf("first parse: %s %s\n",row->term[0],row->term[1]);
      struct split_text* bp1_info = split_line(row->term[0],":");
      struct split_text* bp2_info = split_line(row->term[1],":");
      if(bp1_info->cols == 3 && bp2_info->cols == 3) {
	struct split_text* range1 = split_line(bp1_info->term[2],"-");
	struct split_text* range2 = split_line(bp2_info->term[2],"-");
	int stem_start1 = atoi(range1->term[0]);
	int stem_end1 = atoi(range1->term[1]);
	int length1 = stem_end1 - stem_start1 + 1;
	int stem_start2 = atoi(range2->term[0]);
	int stem_end2 = atoi(range2->term[1]);
	int length2 = stem_end2 - stem_start2 + 1;
	if(length1 == length2) {
	  size += length1;
	} else {
	  printf("lengths don't match in stem file.\n");
	  exit(1);
	}
      }
    } else {
      printf("unexpected number of columns.\n");
      exit(1);
    }
    free_split_text(row);
  }
  return size;
}

struct bp_array* read_basepair_file(char *basepair_file) {
  char line[MAX_LINE_LENGTH];

  FILE* F = fopen(basepair_file,"r");
  if(F == NULL) {
    printf("Error opening the file:\"%s\".\n Check permissions/file name.\n",basepair_file);
    exit(1);
  }

  printf("reading basepair file: %s\n",basepair_file);
  int size = count_basepairs(F);
  printf("allocating bp_array of size: %d\n",size);
  struct bp_array* bp_array = create_bp_array(size);
  fseek(F,0,0); // reset file

  FILE* OUT = fopen("pl_basepairs.txt","w");
  int count = 0;
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    struct split_text* row = split_line(line," ");
    if(row->cols >= 4) {
      if((strcmp(row->term[2],"W")==0)&&(strcmp(row->term[3],"W")==0)&&(strcmp(row->term[4],"C")==0)) {
	//printf("first parse: %s %s\n",row->term[0],row->term[1]);
	struct split_text* bp1_info = split_line(row->term[0],":");
	struct split_text* bp2_info = split_line(row->term[1],":");
	if(bp1_info->cols == 2 && bp2_info->cols == 2) {
	  struct bp* bp = create_bp(bp1_info->term[0][0],bp2_info->term[0][0],atoi(bp1_info->term[1]),atoi(bp2_info->term[1]));
	  set_bp(bp_array,bp,count);
	} else if(bp1_info->cols == 3 && bp2_info->cols == 3) {
	  // access first character/convert to char with ->term[0][0]
	  fprintf(OUT,"%c %d %c %d\n",bp1_info->term[0][0],atoi(bp1_info->term[2]),bp2_info->term[0][0],atoi(bp2_info->term[2]));
	  struct bp* bp = create_bp(bp1_info->term[0][0],bp2_info->term[0][0],atoi(bp1_info->term[2]),atoi(bp2_info->term[2]));
	  set_bp(bp_array,bp,count);	  
	} else {
	  printf("Unexpected situation in basepair file: %s",basepair_file);
	}
	free_split_text(bp1_info);
	free_split_text(bp2_info);
	count++;
      }
    }
    free_split_text(row);
  }
  printf("reading basepair file complete.\n");
  fclose(OUT);
  return bp_array;
}

struct bp_array* read_stems_file(char *stems_file) {
  char line[MAX_LINE_LENGTH];

  FILE* F = fopen(stems_file,"r");
  if(F == NULL) {
    printf("Error opening the file:\"%s\".\n Check permissions/file name.\n",stems_file);
    exit(1);
  }

  printf("reading stems file: %s\n",stems_file);
  int size = count_stem_basepairs(F);
  printf("allocating bp_array of size: %d\n",size);
  struct bp_array* bp_array = create_bp_array(size);
  fseek(F,0,0); // reset file

  FILE* OUT = fopen("pl_basepairs.txt","w");
  int count = 0;
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    struct split_text* row = split_line(line," ");
    if(row->cols == 2) {
      //printf("first parse: %s %s\n",row->term[0],row->term[1]);
      struct split_text* stem_info1 = split_line(row->term[0],":");
      struct split_text* stem_info2 = split_line(row->term[1],":");
      if(stem_info1->cols == 3 && stem_info2->cols == 3) {
	struct split_text* range1 = split_line(stem_info1->term[2],"-");
	struct split_text* range2 = split_line(stem_info2->term[2],"-");
	int stem_start1 = atoi(range1->term[0]);
	int stem_end1 = atoi(range1->term[1]);
	int stem_start2 = atoi(range2->term[0]);
	int stem_end2 = atoi(range2->term[1]);
	for(int pos1 = stem_start1; pos1 <= stem_end1; pos1++) {
	  int pos2 = stem_end2 - (pos1 - stem_start1);
	  // for chain_id, convert from string to char with ->term[0][0] i.e. access first char
	  fprintf(OUT,"%c %d %c %d\n",stem_info1->term[0][0],pos1,stem_info2->term[0][0],pos2);
	  struct bp* bp = create_bp(stem_info1->term[0][0],stem_info2->term[0][0],pos1,pos2);
	  set_bp(bp_array,bp,count);
	  count++;
	}
      } else {
	printf("Unexpected situation in basepair file: %s",stems_file);
      }
      free_split_text(stem_info1);
      free_split_text(stem_info2);
    } else {
      printf("error");
      exit(1);
    }
    free_split_text(row);
  }
  printf("reading basepair file complete.\n");
  fclose(OUT);
  return bp_array;
}

/*
 * MOTIF INPUT
 */

struct motif_group* read_motif_file(char* motif_file) {
  char line[MAX_LINE_LENGTH];

  FILE* F = fopen(motif_file,"r");
  if(F == NULL) {
    printf("Error opening the file:\"%s\".\n Check permissions/file name.\n",motif_file);
    exit(1);
  }  

  struct motif_group* motif_group = create_motif_group(1);
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    struct split_text* row = split_line(line,"-");
    struct motif* motif = create_motif(row->cols);
    for(int i=0; i<row->cols; i++) {
      set_motif_loop(motif,i,atoi(row->term[i]));
    }
    append_motif(motif_group,motif);
  }
  fclose(F);
  return motif_group;
}
