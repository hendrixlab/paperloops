#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "RNA_molecule.h"
#include "RNA_PDB.h"

#define MAX_LINE_LENGTH 1000

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
  return nuc_name;
}

int get_model_count(FILE* F) {
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

int* get_model_locs(FILE* F, int* model_count, int** model_locs) {
  char line[MAX_LINE_LENGTH];
  char rec_type[6]; // the max length for a record type is 6
  int FOUND_MODEL = 0;
  int FOUND_ATOM = 0;

  // get the number of models 
  *model_count = get_model_count(F);

  // reset file and allocate memory
  fseek(F,0,0);
  *model_locs = malloc(*model_count * sizeof(int));

  // read the model locations
  while(fgets(line,MAX_LINE_LENGTH,F)) {
    strncpy(rec_type,line,6);
    if(strcmp(rec_type,"ATOM  ") == 0) {
      if(!FOUND_MODEL && !FOUND_ATOM) {
	FOUND_ATOM = 1;
	int model_num = 1;
	int offset = ftell(F) - strlen(line);
	printf("model %d offset: %d\n",model_num,offset);
	(*model_locs)[model_num-1] = offset;	
      }
    } else if(strcmp(rec_type,"MODEL ") == 0) {
      // MODEL        1                                                                  
      FOUND_MODEL = 1;
      char model_num_str[5];
      // model serial number at positions 11-14
      sprintf(model_num_str,"%.4s",line+11);
      int model_num = atoi(model_num_str);
      int offset = ftell(F);
      printf("model %d offset: %d\n",model_num,offset);
      (*model_locs)[model_num-1] = offset;
    }
  }
  //*model_locs = temp_model_locs;
}

struct rna** read_pdb_file(char* file_name, int* model_count) {
  FILE* F;
  char line[MAX_LINE_LENGTH];
  char rec_type[6]; // the max length for a record type is 6
  
  F = fopen(file_name,"r");
  if(F == NULL) {
    printf("Error opening the file. Check permissions/file name.\n");
    exit(1);
  }

  // get the number and locations of each model
  int* model_locs;
  get_model_locs(F, model_count, &model_locs);

  struct rna** models = malloc(sizeof(struct rna*) * (*model_count));
  // loop through the model locations and read in data: 
  for(int i=0; i< *model_count; i++) {
    fseek(F,model_locs[i],0);  // jump to the first ATOM line of model
    struct rna* RNA = create_rna();
    while(fgets(line,MAX_LINE_LENGTH,F)) {
      strncpy(rec_type,line,6);
      if(strcmp(rec_type,"ATOM  ") == 0) {	
	int nuc_num = extract_nuc_num(line);
	char* nuc_name = extract_nuc_name(line);	
	if(!last_nuc_has_num(RNA,nuc_num)) {
	  struct nuc* n = create_nuc(nuc_num,nuc_name);
	  add_nuc(RNA,n);
	}
	int atom_num = extract_atom_num(line);
	char* atom_name = extract_atom_name(line);
	double* r = extract_coord(line);
	struct atom* a = create_atom(atom_num,atom_name,r);
	add_atom(RNA,a);
      }
    }
    models[i] = RNA;
  }
  fclose(F);
  return models;
}
