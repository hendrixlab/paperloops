/*
 * This library contains the definitions of structures and functions reading PDB files
 * for RNA structures and storing into the RNA data structure defined in RNA_molecule.c
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "motifs.h"
#include "rna_io.h"
#include "rna_data.h"
#include "transform.h"

#define MAX_LINE_LENGTH 1000


/* FUNCTIONS */

int is_numeric(char* string) {
  char* ptr = string;
  while(*ptr != '\0') {
    if(!isdigit((unsigned char)*ptr)) {
      return 0;
      break;
    }
    ptr++;
  }
  return 1;
}

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
  // should be pos 18-20 in 1-based (3 char + \0 gives 4)
  char* nuc_name = malloc(4 * sizeof(char));
  
  int count = 0;
  for(int i=0; i<3; i++) {
    char c = line[17 + i];
    if(c != ' ' && c != '\0') {      
      nuc_name[count++] = c;
    }
  }
  nuc_name[count] = '\0';
  
  if(count == 0) {
    strcpy(nuc_name,"N");
  }
  
  // need to map PSU to U and other HETATM nucleotides
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
      //if((strncmp(line,"ATOM  ",6) == 0)||(strncmp(line,"HETATM",6) == 0)) {
      if(strncmp(line,"ATOM  ",6) == 0) {
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
          // create a nucleotide
          //printf("adding %s to %d\n",nuc_name,chain_length(get_last_chain(RNA_model)));
	  struct nuc* n = create_nuc(nuc_num,nuc_name);
	  add_nuc(RNA_model,n);	  
	  //printf("Adding: %c %d %s\n",chain_id,nuc_num,nuc_name);
	  //printf("Current: %s\n",get_chain_sequence(get_last_chain(RNA_model)));
	} else {
          free(nuc_name);
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
    struct chain* chain = get_first_chain(RNA_model);
    int total_chains = 0;
    while(chain) {
      set_chain_sequence(chain);
      total_chains++;
      chain = next_chain(chain);      
    }
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
        free_split_text(range1);
        free_split_text(range2);
	if(length1 == length2) {
	  size += length1;
	} else {
	  printf("lengths don't match in stem file.\n");
	  exit(1);
	}
      } else if(bp1_info->cols == 2 && bp2_info->cols == 2) {
	struct split_text* range1 = split_line(bp1_info->term[1],"-");
	struct split_text* range2 = split_line(bp2_info->term[1],"-");
	int stem_start1 = atoi(range1->term[0]);
	int stem_end1 = atoi(range1->term[1]);
	int length1 = stem_end1 - stem_start1 + 1;
	int stem_start2 = atoi(range2->term[0]);
	int stem_end2 = atoi(range2->term[1]);
	int length2 = stem_end2 - stem_start2 + 1;
        free_split_text(range1);
        free_split_text(range2);
	if(length1 == length2) {
	  size += length1;
	} else {
	  printf("lengths don't match in stem file.\n");
	  exit(1);
	}
      } else {
	printf("parsing error in stems file.\n");
	exit(1);
      }
      free_split_text(bp1_info);
      free_split_text(bp2_info);
    } else {
      printf("unexpected number of columns in stems file.\n");
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

  //FILE* OUT = fopen("pl_basepairs.txt","w");
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
	  set_bp(bp_array,count,bp);
	} else if(bp1_info->cols == 3 && bp2_info->cols == 3) {
	  // access first character/convert to char with ->term[0][0]
	  //fprintf(OUT,"%c %d %c %d\n",bp1_info->term[0][0],atoi(bp1_info->term[2]),bp2_info->term[0][0],atoi(bp2_info->term[2]));
	  struct bp* bp = create_bp(bp1_info->term[0][0],bp2_info->term[0][0],atoi(bp1_info->term[2]),atoi(bp2_info->term[2]));
	  set_bp(bp_array,count,bp);
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
  //fclose(OUT);
  fclose(F);
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
        free_split_text(range1);
        free_split_text(range2);
	for(int pos1 = stem_start1; pos1 <= stem_end1; pos1++) {
	  int pos2 = stem_end2 - (pos1 - stem_start1);
	  // for chain_id, convert from string to char with ->term[0][0] i.e. access first char
	  char chain_id1 = stem_info1->term[0][0];
	  char chain_id2 = stem_info2->term[0][0];
	  fprintf(OUT,"%c %d %c %d\n",chain_id1,pos1,chain_id2,pos2);
	  struct bp* bp = create_bp(chain_id1,chain_id2,pos1,pos2);
	  set_bp(bp_array,count,bp);
	  count++;
	}
      } else if(stem_info1->cols == 2 && stem_info2->cols == 2) {
	struct split_text* range1 = split_line(stem_info1->term[1],"-");
	struct split_text* range2 = split_line(stem_info2->term[1],"-");
	int stem_start1 = atoi(range1->term[0]);
	int stem_end1 = atoi(range1->term[1]);
	int stem_start2 = atoi(range2->term[0]);
	int stem_end2 = atoi(range2->term[1]);
        free_split_text(range1);
        free_split_text(range2);
	for(int pos1 = stem_start1; pos1 <= stem_end1; pos1++) {
	  int pos2 = stem_end2 - (pos1 - stem_start1);
	  // for chain_id, convert from string to char with ->term[0][0] i.e. access first char
	  char chain_id1 = stem_info1->term[0][0];
	  char chain_id2 = stem_info2->term[0][0];
	  fprintf(OUT,"%c %d %c %d\n",chain_id1,pos1,chain_id2,pos2);
	  struct bp* bp = create_bp(chain_id1,chain_id2,pos1,pos2);
	  set_bp(bp_array,count,bp);
	  count++;
	}
      } else {
	printf("Unexpected situation in stems file: %s",stems_file);
	exit(1);
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
  fclose(F);
  return bp_array;
}

/*
 * MOTIF INPUT
 */


int read_max_motif_length(FILE* MOTIF_FILE) {
  char line[MAX_LINE_LENGTH];

  int L_max = 0;
  while(fgets(line,MAX_LINE_LENGTH,MOTIF_FILE)) {
    struct split_text* row = split_line(line,"-");
    int n = row->cols;   // row->cols is n 
    struct motif* motif = create_motif(n);
    for(int i=0; i<row->cols; i++) {  // motif has to be built with loop lengths
      set_motif_loop(motif,i,atoi(row->term[i]));
    }
    int L = get_motif_sum(motif);
    if(L > L_max) {
      L_max = L;
    }
    free_motif(motif);
  }
  return L_max;
}

struct motif* read_input_motif(char* input_motif) {
  struct motif* motif;
  if(input_motif[0] == '-') {
    motif = create_motif(1);
    char* substring = input_motif + 1;
    if(is_numeric(substring)) {
      set_motif_loop(motif,0,atoi(substring));
    } else {
      printf("Error: incorrect motif formatting: %s\n", input_motif);
      exit(0);
    }
  } else {
    struct split_text* row = split_line(input_motif,"-");
    int n = row->cols;
    if(n == 1) {
      // if only one number, then n=0 (no helices)
      motif = create_motif(0);
      if(is_numeric(row->term[0])) {
        set_motif_loop(motif,0,atoi(row->term[0]));
      } else {
        printf("Error: incorrect motif formatting: %s\n", input_motif);
        exit(0);
      }
    } else {
      // should be any number of loops/helices n
      motif = create_motif(n);
      for(int i=0; i<row->cols; i++) {  // motif has to be built with loop lengths
        if(is_numeric(row->term[0])) {
          set_motif_loop(motif,i,atoi(row->term[i]));
        } else {
          printf("Error: incorrect motif formatting: %s\n", input_motif);
          exit(0);
        }
      }
    }
  }
  int L = get_motif_sum(motif);
  struct hash_table* motif_hash = create_hash_table(L); // not sure what size to use here
  long int id = compute_motif_id(motif,L);
  printf("staring id: %d\n",id);
  set_motif_id(motif,id);
  struct motif* curr_motif = motif;
  hash_table_set(motif_hash,id,curr_motif); // each rotation and starting motif gets added to motif_hash

  for(int j=0; j<motif->num; j++) { // apart from repeats, motif with n loops has at most n rotations
    struct motif* prev_motif = curr_motif; // initially prev_motif has motif that was read in, but is updated
    curr_motif = rotate_motif(curr_motif,L); 
    long int id = get_motif_id(curr_motif);
    printf("rotated id: %d\n",id);

    if(id < 0) {
      printf("ERROR in motif ID: %l\t%d\n",id,curr_motif->num);
      exit(0);
    }
    if(hash_table_empty(motif_hash,id)) { // only add if it is not visited/stored yet
      hash_table_set(motif_hash,id,curr_motif); // each rotation gets added to motif_hash
      set_next_motif(prev_motif,curr_motif); // attach the rotation so motif_lib has all of them
      printf("Attaching motif rotation %d: ",j);
      print_motif(curr_motif,'\n');
    }
  }
  return motif;
}

struct motif_group* read_motif_file(char* motif_file) {
  char line[MAX_LINE_LENGTH];

  FILE* MOTIF_FILE = fopen(motif_file,"r");
  if(MOTIF_FILE == NULL) {
    printf("Error opening the file:\"%s\".\n Check permissions/file name.\n",motif_file);
    exit(1);
  }  

  struct motif_group* motif_group = create_motif_group(1);
  int L_max = read_max_motif_length(MOTIF_FILE);
  printf("Read in max motif length: %d\n",L_max);
  fseek(MOTIF_FILE,0,0); // reset file
  struct hash_table* motif_hash = create_hash_table(L_max); // not sure what size to use here
  
  while(fgets(line,MAX_LINE_LENGTH,MOTIF_FILE)) {
    if(line[0] != '#') {
      struct motif* motif;
      if(line[0] == '-') {
        motif = create_motif(1);
        char* substring = line + 1;
        if(is_numeric(substring)) {
          set_motif_loop(motif,0,atoi(substring));
        } else {
          printf("Error: incorrect motif formatting: %s\n", line);
          exit(0);
        }
      } else {
        struct split_text* row = split_line(line,"-");
        int n = row->cols;
        if(n == 1) {
          // if only one number, then n=0 (no helices)
          motif = create_motif(0);
          if(is_numeric(row->term[0])) {
            set_motif_loop(motif,0,atoi(row->term[0]));
          } else {
            printf("Error: incorrect motif formatting: %s\n", line);
            exit(0);
          }
        } else {
          // should be any number of loops/helices n
          motif = create_motif(n);
          for(int i=0; i<row->cols; i++) {  // motif has to be built with loop lengths
            if(is_numeric(row->term[0])) {
              set_motif_loop(motif,i,atoi(row->term[i]));
            } else {
              printf("Error: incorrect motif formatting: %s\n", line);
              exit(0);
            }
          }
        }
      }
      long int id = compute_motif_id(motif,L_max);
      set_motif_id(motif,id);
      hash_table_set(motif_hash,id,motif);
      append_motif(motif_group,motif);
      // motif added, now add all the rotations
      for(int j=0; j<motif->num; j++) { // apart from repeats, motif with n loops has at most n rotations
        struct motif* prev_motif = motif; // initially prev_motif has motif that was read in, but is updated
        motif = rotate_motif(motif,L_max); 
        long int id = get_motif_id(motif);
        if(id < 0) {
          printf("ERROR in motif ID: %l\t%d\n",id,motif->num);
          exit(0);
        }
        if(hash_table_empty(motif_hash,id)) { // only add if it is not visited/stored yet
          hash_table_set(motif_hash,id,motif); // each rotation gets added to motif_hash
          set_next_motif(prev_motif,motif); // attach the rotation so motif_lib has all of them
        }
      }
    }
  }
  free_hash_table(motif_hash);
  fclose(MOTIF_FILE);
  return motif_group;
}


/***********************************************************
 * processing input biophysical data to compute transforms *
 ***********************************************************/

int valid_dinucleotide(struct nuc* n1, struct nuc* n2) {
  struct atom* a1 = get_last_backbone_atom(n1);
  struct atom* a2 = get_first_backbone_atom(n2);
  if(a1 && a2) {
    if((strcmp(atom_name(a1),"O3'") == 0)&&(strcmp(atom_name(a2),"P") == 0)) {
      double* diff = subtract(atom_coord(a1),atom_coord(a2));
      double d = norm(diff);
      free_vector(diff);
      if(d < 2) {
	return 1;
      }    
    } else {
      printf("Unexpected atoms found for nucs %d and %d: %s %s\n",nuc_num(n1),nuc_num(n2),atom_name(a1),atom_name(a2));
      exit(1);
    }
  }
  return 0;
}

struct transform_library* build_backbone_transform_library(struct rna_data* RNA) {
  assert(RNA);
  if(get_model_count(RNA) == 0) {
    printf("No models found for RNA.\n");
    exit(1);
  }
  double d_sum = 0.0;
  int d_tot = 0;
  struct transform_library* backbone_library = create_transform_library();
  for(int i=0; i<get_model_count(RNA); i++) {
    struct rna_model* RNA_model = get_rna_model(RNA,i); 
    struct chain* chain = get_first_chain(RNA_model);
    while(chain) {
      struct nuc* n1 = get_first_nuc(chain);
      struct nuc* n2 = next_nuc(chain,n1);
      while(before_end(chain,n2)) { // while n2 is not a sentinel & is a valid nt
	if(valid_dinucleotide(n1,n2)) { // checks that no gap b/w n1 and n2
	  struct frame* f1 = compute_frame(n1);
	  struct frame* f2 = compute_frame(n2);
	  struct transform* T = compute_transform(f1,f2);
	  // for calculating the avg displacement
	  double d = get_transform_displacement(T);
	  d_sum += d;
	  d_tot++;
	  // store backbone in transform_library
	  append_transform(backbone_library,T);
	  free_frame(f1); // free memory of each transform after use
	  free_frame(f2);
	}
	n1 = n2;
	n2 = next_nuc(chain,n1);	   
      }
      chain = next_chain(chain);
    }
  }
  printf("Found %d dinucleotide pairs.\n",transform_library_size(backbone_library));
  printf("Avg backbone displacement: %.2f\n",d_sum/d_tot);
  return backbone_library;
}

struct transform_library* build_basepair_transform_library(struct rna_data* RNA, struct bp_array* bp_array) {
  assert(RNA);
  if(get_model_count(RNA) == 0) {
    printf("No models found for RNA.\n");
    exit(1);
  }
  double d_sum = 0.0;
  int d_tot = 0;
  struct transform_library* basepair_library = create_transform_library();
  for(int i=0; i<get_model_count(RNA); i++) {      
    struct rna_model* RNA_model = get_rna_model(RNA,i);
    for(int j=0; j<get_bp_array_size(bp_array);j++) {
      struct bp* bp = get_bp(bp_array,j);
      char chain_id1 = get_bp_chain_id(bp,1);
      char chain_id2 = get_bp_chain_id(bp,2);
      //printf("checking chain: %c %c\n",chain_id1,chain_id2);
      struct chain* chain1 = get_chain_by_id(RNA_model,chain_id1);
      struct chain* chain2 = get_chain_by_id(RNA_model,chain_id2);
      int pos1 = get_bp_pos(bp,1);
      int pos2 = get_bp_pos(bp,2);
      //printf("checking pos: %d %d\n",pos1,pos2);
      struct nuc* n1 = get_nuc_by_idx(chain1,pos1);
      struct nuc* n2 = get_nuc_by_idx(chain2,pos2);
      //printf("basepair: %d %d\n",nuc_num(n1),nuc_num(n2));
      if(n1 && n2) {
	struct frame* f1 = compute_frame(n1);
	struct frame* f2 = compute_frame(n2);
	struct transform* T = compute_transform(f1,f2);
	// for calculating the avg displacement
	double d = get_transform_displacement(T);
	d_sum += d;
	d_tot++;      
	append_transform(basepair_library,T);
	free_frame(f1);
	free_frame(f2);
      }
    }
  }
  printf("found %d base pairs.\n",transform_library_size(basepair_library));
  printf("Avg basepair displacement: %.2f\n",d_sum/d_tot);
  return basepair_library;
}

struct library* build_library(char* pdb_file, char* stems_file) {
  // build combined library of bp and bb transforms
  
  printf("reading in PDB file: %s\n",pdb_file);
  struct rna_data* RNA = read_pdb_file(pdb_file);

  // building libraries
  struct library* rna_lib = create_library();
  printf("Building backbone libraries:\n");
  rna_lib->bb = build_backbone_transform_library(RNA);
  
  printf("reading in stems file: %s\n",stems_file);
  struct bp_array* basepairs = read_stems_file(stems_file);

  printf("Building basepair libraries:\n");  
  rna_lib->bp = build_basepair_transform_library(RNA,basepairs);

  // memory clean-up
  free_bp_array(basepairs); // maybe want to save for later? 
  free_rna_data(RNA);
  
  return rna_lib;
}

/******************
 * HELIX FUNCTION *
 ******************/

struct library* build_helix_library(char* helix_pdb_file) {
  // build combined library of bp and bb transforms

  struct library* helix_lib = create_library();

  if(helix_pdb_file) {
    struct rna_data* helix = read_pdb_file(helix_pdb_file);

    // stems file:
    char helix_stems_file[100];
    strcpy(helix_stems_file,helix_pdb_file);
    strcat(helix_stems_file,".stems.txt");
    printf("reading in helix stems file: %s\n",helix_stems_file);
    struct bp_array* helix_basepairs = read_stems_file(helix_stems_file);
    
    helix_lib->bb = build_backbone_transform_library(helix);
    helix_lib->bp = build_basepair_transform_library(helix,helix_basepairs);

    free_bp_array(helix_basepairs);
    free_rna_data(helix);
  }
  return helix_lib;
}
