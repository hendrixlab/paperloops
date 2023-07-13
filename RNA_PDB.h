/*
 * This library contains the definitions of structures and functions reading PDB files
 * for RNA structures and storing into the RNA data structure defined in RNA_molecule.c
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#ifndef _RNA_PDB
#define _RNA_PDB

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "RNA_PDB.h"
#include "RNA_molecule.h"

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

#endif
