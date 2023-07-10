#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "RNA_molecule.h"

#define MAX_LINE_LENGTH 1000

/* FUNCTION PROTOTYPES */

double* extract_coord(char* line);

double* extract_coord(char* line);

int extract_atom_num(char* line);

char* extract_atom_name(char* line);

int extract_nuc_num(char* line);

char* extract_nuc_name(char* line);

int get_model_count(FILE* F);

int* get_model_locs(FILE* F, int* model_count, int** model_locs);

struct rna** read_pdb_file(char* file_name, int* model_count);

