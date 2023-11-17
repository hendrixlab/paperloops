/*
 * This library contains the definitions of structures and functions reading PDB files
 * for RNA structures and storing into the RNA data structure defined in RNA_molecule.c
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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
  char rec_type[6]; // the max length for a record type is 6
  
  int model_count = 1; // one model by default 
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    strncpy(rec_type,line,6);
    if(strcmp(rec_type,"MODEL ") == 0) {
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
  char rec_type[6]; // the max length for a record type is 6
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
    strncpy(rec_type,line,6);
    if(strcmp(rec_type,"ATOM  ") == 0) {
      if(!FOUND_MODEL && !FOUND_ATOM) {
	FOUND_ATOM = 1;
	int model_num = 0;
	int offset = ftell(F) - strlen(line);
	set_model_loc(RNA,model_num,offset);
      }
    } else if(strcmp(rec_type,"MODEL ") == 0) {
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
  char rec_type[6]; // the max length for a record type is 6
  fseek(F,get_model_loc(RNA,model_num),0);  // jump to the first ATOM line of model
  char prev_chain_id = '?';
  int COUNT1 = 0;
  int COUNT2 = 0; // second method of counting: double checking
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    strncpy(rec_type,line,6);
    if(strcmp(rec_type,"ATOM  ") == 0) {
      char chain_id = extract_chain_id(line);
      if(prev_chain_id == '?') {
	prev_chain_id = chain_id;
	COUNT1++; // add to count for first chain
      }
      if(prev_chain_id != chain_id) {
	COUNT1++; // add for all other chains
	prev_chain_id = chain_id;
      }
    } else if(strcmp(rec_type,"TER   ")) {
      COUNT2++;
    } else if(strcmp(rec_type,"ENDMDL") == 0) {
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
  char rec_type[6]; // the max length for a record type is 6
  
  FILE* F = fopen(input_file_name,"r");
  if(F == NULL) {
    printf("Error opening the file:\"%s\".\n Check permissions/file name.\n",input_file_name);
    exit(1);
  }

  struct rna_data* RNA = create_rna_data(input_file_name);
  // get the number and locations of each model, allocate memory
  read_model_locs(F, RNA);
  allocate_rna_models(RNA);
  
  for(int model_num=0; model_num< get_model_count(RNA); model_num++) {
    //int chain_count = read_model_chain_count(FILE,RNA,model_num); // get the number of chains for this model
    struct rna_model* RNA_model = create_model(); // create RNA_model
    // loop through the model locations and read in data: 
    fseek(F,get_model_loc(RNA,model_num),0);  // jump to the first ATOM line of model
    while(fgets(line,MAX_LINE_LENGTH,F)) {
      strncpy(rec_type,line,6);
      if((strcmp(rec_type,"ATOM  ") == 0)||(strcmp(rec_type,"HETATM") == 0)) {
	//if(strcmp(rec_type,"ATOM  ") == 0) {
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
      } else if(strcmp(rec_type,"TER   ") == 0) {
	// expecting end of chain
	;
      } else if(strcmp(rec_type,"ENDMDL") == 0) {
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

struct split_text*  split_line(char *line, char *delim) {
  struct split_text* row = malloc(sizeof(struct split_text));
  row->term = NULL;
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

  int count = 0;
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    struct split_text* row = split_line(line," ");
    if(row->cols >= 4) {
      if((strcmp(row->term[2],"W")==0)&&(strcmp(row->term[3],"W")==0)&&(strcmp(row->term[4],"C")==0)) {
	//printf("first parse: %s %s\n",row->term[0],row->term[1]);
	struct split_text* bp1_info = split_line(row->term[0],":");
	struct split_text* bp2_info = split_line(row->term[1],":");
	//printf("%s %s %s %s\n",bp1_info->term[0],bp2_info->term[0],bp1_info->term[1],bp2_info->term[1]);
	if(bp1_info->cols == 2 && bp2_info->cols == 2) {
	  struct bp* bp = create_bp(bp1_info->term[0][0],bp2_info->term[0][0],atoi(bp1_info->term[1]),atoi(bp2_info->term[1]));
	  set_bp(bp_array,bp,count);
	} else if(bp1_info->cols == 3 && bp2_info->cols == 3) {
	  struct bp* bp = create_bp(bp1_info->term[0][0],bp2_info->term[0][0],atoi(bp1_info->term[2]),atoi(bp2_info->term[2]));
	  set_bp(bp_array,bp,count);	  
	} else {
	  printf("Unexpected situation in basepair file: %s",basepair_file);
	}
	count++;
      }
    }
  }
  printf("reading basepair file complete.\n");
  return bp_array;
}
