/*
 * This library contains the definitions of structures and functions for implementing
 * the RNA data structure
 *
 * author: David Hendrix david.hendrix@oregonstate.edu
 */

#ifndef _RNA_DATA
#define _RNA_DATA

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "vectorfun.h"
#include "rna_data.h"

/* Primary RNA data structure */
struct rna_data;

/* Structure used to represent RNA molecule linked list */
struct rna_model;

/* Structure used to represent the nucleotide node */
struct nuc;

/* Structure used to represent the atom node */
struct atom;

struct rna_data* create_rna_data(char* input_file_name);
struct rna_model* create_rna_model();
struct nuc* create_nuc(int nuc_num, char* nuc_name);
struct atom* create_atom(int atom_num, char* atom_name, double* r);
void free_rna_model(struct rna_model* RNA_model);
void free_nuc(struct nuc* n);
void free_atom(struct atom* a);
int rna_model_is_empty(struct rna_model* RNA_model);
int nuc_backbone_is_empty(struct nuc* n);
int nuc_sugar_is_empty(struct nuc* n);
int nuc_base_is_empty(struct nuc* n);
void set_model_count(struct rna_data* RNA, int model_count);
int get_model_count(struct rna_data* RNA);
void allocate_model_locs(struct rna_data* RNA);
void set_model_loc(struct rna_data* RNA, int model_num, int offset);
int get_model_loc(struct rna_data* RNA, int model_num);
void allocate_rna_models(struct rna_data* RNA);
void set_rna_model(struct rna_data* RNA, int model_num, struct rna_model* RNA_model);
void add_nuc(struct rna_model* RNA_model, struct nuc* n);
void add_atom(struct rna_model* RNA_model, struct atom* a);
int is_backbone_atom(struct atom* a);
int is_sugar_atom(struct atom* a);
void add_backbone_atom(struct rna_model* RNA_model, struct atom* a);
void add_sugar_atom(struct rna_model* RNA_model, struct atom* a);
void add_base_atom(struct rna_model* RNA_model, struct atom* a);
int new_nuc_num(struct rna_model* RNA_model, int nuc_num);
struct nuc* get_first_nuc(struct rna_model* RNA_model);
struct nuc* get_nuc_by_num(struct rna_model* RNA_model, int nuc_num);
int before_end(struct rna_model* RNA_model, struct nuc* n);
struct nuc* next_nuc(struct rna_model* RNA_model, struct nuc* n);
struct nuc* prev_nuc(struct rna_model* RNA_model, struct nuc* n);
struct atom* prev_nuc_last_backbone_atom(struct rna_model* RNA_model, struct nuc* n);
struct atom* next_nuc_first_backbone_atom(struct rna_model* RNA_model, struct nuc* n);
struct atom* get_first_backbone_atom(struct nuc* n);
struct atom* get_first_sugar_atom(struct nuc* n);
struct atom* next_atom(struct atom* a);
struct atom* get_atom_by_num(struct nuc* n, int atom_num);
struct atom* print_sugar_atoms(struct nuc* n);
struct atom* get_sugar_atom_by_name(struct nuc* n, char* atom_name);
struct atom* remove_sugar_atom_by_name(struct nuc* n, char* atom_name);
int get_nuc_num(struct nuc* n);
char* get_nuc_name(struct nuc* n);
int get_atom_num(struct atom* a);
char* get_atom_name(struct atom* a);
char* get_atom_type(struct atom* a);
double torsion(struct atom* p1, struct atom* p2, struct atom* p3, struct atom *p4);
void compute_backbone_torsions(struct rna_data* RNA_data);
void compute_sugar_torsions(struct rna_data* RNA_data);
void compute_model_backbone_torsions(struct rna_model* RNA_model);
void compute_model_sugar_torsions(struct rna_model* RNA_model);
void print_backbone_torsions(struct rna_data* RNA);
void print_sugar_torsions(struct rna_data* RNA);
char* set_rna_model_sequence(struct rna_model* RNA_model);
char* get_rna_model_sequence(struct rna_model* RNA_model);

#endif
