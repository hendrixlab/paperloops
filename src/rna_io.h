/*
 * This library contains the definitions of structures and functions reading PDB files
 * for RNA structures and storing into the RNA data structure defined in RNA_molecule.c
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 */

#ifndef _RNA_IO
#define _RNA_IO

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "rna_io.h"
#include "rna_data.h"

/*******************
 * data structures *
 *******************/

struct split_text {
  char** term;
  int cols;
};

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

struct split_text* create_split_text();

void free_split_text(struct split_text* row);

int get_split_text_cols(struct split_text* row);

char* get_split_text_term(struct split_text* row, int i);

struct split_text* split_line(char *line, char *delim);

int count_basepairs(FILE* F);

struct bp_array* read_basepair_file(char *stems_file);

struct bp_array* read_stems_file(char *stems_file);

struct motif* read_input_motif(char* input_motif);

struct motif_group* read_motif_file(char* motif_file);

/***********************************************************
 * processing input biophysical data to compute transforms *
 ***********************************************************/

int valid_dinucleotide(struct nuc* n1, struct nuc* n2);

struct transform_library* build_backbone_transform_library(struct rna_data* RNA);

struct transform_library* build_basepair_transform_library(struct rna_data* RNA, struct bp_array* bp_array);

struct library* build_library(char* pdb_file, char* stems_file);

struct library* build_helix_library(char* helix_pdb_file);

#endif
