/*
 * This library contains the definitions of structures and functions reading PDB files
 * for RNA structures and storing into the RNA data structure defined in RNA_molecule.c
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#ifndef _RNA_IO
#define _RNA_IO

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "rna_io.h"
#include "rna_data.h"

/* DATA TYPES */

struct split_text;

/* FUNCTION PROTOTYPES */

double* extract_coord(char* line);

double* extract_coord(char* line);

int extract_atom_num(char* line);

char* extract_atom_name(char* line);

int extract_nuc_num(char* line);

char* extract_nuc_name(char* line);

char extract_chain_id(char* line);

int read_model_count(FILE* F);

int* read_model_locs(FILE* F, struct rna_data* RNA);

struct rna_data* read_pdb_file(char* file_name);

struct split_text*  split_line(char *line, char *delim);

int count_basepairs(FILE* F);

struct bp_array* read_basepair_file(char *stems_file);

#endif
