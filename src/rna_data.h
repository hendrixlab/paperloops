/*
 * This library contains the definitions of structures and functions for implementing
 * the RNA data structure
 *
 * (C) 2023 David Hendrix david.hendrix@oregonstate.edu
 */

#ifndef _RNA_DATA
#define _RNA_DATA

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "vectorfun.h"
#include "rna_data.h"


/*
 *  Data Structure Definitions 
 *  includes nodes for linked lists describing:
 *     -sequence of nucleotides (linked list of "nuc" nodes   
 *     -collection of atoms. There are three linked lists of atoms:
 *         -backbone atoms
 *         -sugar atoms
 *         -base atoms
 */

struct rna_data {
  struct rna_model** models;
  int model_count;
  int *model_locs;
  char* input_file_name;
};
  
struct rna_model {
  // doubly linked list of RNA polymer chains
  struct chain* first_chain;
  struct chain* last_chain;
  int chain_count;  
};

struct chain {
  struct nuc *sentinel; // circular doubly-linked list of nuclotides
  int length;
  char* seq;
  struct chain *next; // next chain in the model
  struct chain *prev; // previous chain in the model
  char chain_id;
  struct nuc **nuc_array; // array of nucleotides
  int nuc_array_size;
};

// nucleotide node
struct nuc {
  struct atom* backbone_first;
  struct atom* backbone_last;
  struct atom* sugar_first;
  struct atom* sugar_last;
  struct atom* base_first;
  struct atom* base_last;
  int num;
  char* name;
  int size;
  double bbta[6]; // backbone torsion angles
  double sta[5];  // sugar torsion angles 
  struct nuc* next; // next nucleotide in the chain
  struct nuc* prev; // prev nucleotide in the chain
};

// atom node
struct atom {
  char* name;
  double* r;
  int num;
  char* type;
  struct atom* next;
  struct atom* prev;
  struct atom* sugar_next; // pointers just for sugar ring
  struct atom* sugar_prev;  
};

// base pair array
struct bp_array {
  struct bp** data;
  int size;
};

// base pair
struct bp {
  char chain_id1;
  char chain_id2;
  int pos1;
  int pos2;
};


struct rna_data* create_rna_data(char* input_file_name);
struct rna_model* create_model();
struct chain* create_chain(char chain_id);
struct nuc* create_nuc(int nuc_num, char* nuc_name);
struct atom* create_atom(int atom_num, char* atom_name, double* r);
void free_atom(struct atom* a);
void free_nuc(struct nuc* n);
void free_chain(struct chain* chain);
void free_rna_data(struct rna_data* RNA_data);
int model_is_empty(struct rna_model* RNA_model);
int chain_is_empty(struct chain* RNA_chain);
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
struct rna_model* get_rna_model(struct rna_data* RNA, int model_num);
void print_rna_model(struct rna_model* RNA_model);
void set_model_chain_count(struct rna_model* RNA_model, int chain_count);
int get_model_chain_count(struct rna_model* RNA_model);
void allocate_model_chains(struct rna_model* RNA_model);
void add_chain(struct rna_model* RNA_model, struct chain* new_chain);
void add_nuc(struct rna_model* RNA_model, struct nuc* n);
void add_atom(struct rna_model* RNA_model, struct atom* a);
int is_backbone_atom(struct atom* a);
int is_sugar_atom(struct atom* a);
void add_backbone_atom(struct rna_model* RNA_model, struct atom* a);
void add_sugar_atom(struct rna_model* RNA_model, struct atom* a);
void add_base_atom(struct rna_model* RNA_model, struct atom* a);
struct chain* get_first_chain(struct rna_model* RNA_model);
struct chain* get_last_chain(struct rna_model* RNA_model);
struct chain* get_chain_by_id(struct rna_model* RNA_model, char chain_id);
struct chain* next_chain(struct chain* RNA_chain);
int new_chain_id(struct rna_model* RNA_model, char chain_id);
int new_nuc_num(struct rna_model* RNA_model, int nuc_num);
int chain_length(struct chain* chain);
struct nuc* get_first_nuc(struct chain* RNA_chain);
struct nuc* get_nuc_by_num(struct chain* RNA_chain, int nuc_num);
struct nuc* get_nuc_by_idx(struct chain* chain, int idx);
int before_end(struct chain* RNA_chain, struct nuc* n);
struct nuc* next_nuc(struct chain* RNA_chain, struct nuc* n);
struct nuc* prev_nuc(struct chain* RNA_chain, struct nuc* n);
void print_chain(struct chain* chain);
struct atom* prev_nuc_last_backbone_atom(struct chain* RNA_chain, struct nuc* n);
struct atom* next_nuc_first_backbone_atom(struct chain* RNA_chain, struct nuc* n);
struct atom* get_first_backbone_atom(struct nuc* n);
struct atom* get_last_backbone_atom(struct nuc* n);
struct atom* get_backbone_atom_by_name(struct nuc* n, char* atom_name);
struct atom* get_first_sugar_atom(struct nuc* n);
struct atom* get_first_base_atom(struct nuc* n);
struct atom* next_atom(struct atom* a);
struct atom* get_atom_by_num(struct nuc* n, int atom_num);
struct atom* get_sugar_atom_by_name(struct nuc* n, char* atom_name);
struct atom* remove_sugar_atom_by_name(struct nuc* n, char* atom_name);
void print_sugar_atoms(struct nuc* n);
int nuc_num(struct nuc* n);
char* nuc_name(struct nuc* n);
int atom_num(struct atom* a);
char* atom_name(struct atom* a);
char* atom_type(struct atom* a);
double* atom_coord(struct atom* a);
double torsion(struct atom* p1, struct atom* p2, struct atom* p3, struct atom *p4);
void compute_backbone_torsions(struct rna_data* RNA_data);
void compute_sugar_torsions(struct rna_data* RNA_data);
void compute_model_backbone_torsions(struct rna_model* RNA_model);
void compute_model_sugar_torsions(struct rna_model* RNA_model);
void print_backbone_torsions(struct rna_data* RNA);
void print_sugar_torsions(struct rna_data* RNA);
void set_rna_nuc_array(struct rna_data* RNA);
void set_chain_nuc_array(struct chain* RNA_chain);

// bp and bp_array tools

struct bp_array* create_bp_array(int size);
void free_bp_array(struct bp_array* bp_array);
struct bp* create_bp(char chain_id1, char chain_id2, int pos1, int pos2);
void free_bp(struct bp* bp);
int get_bp_array_size(struct bp_array* bp_array);
void set_bp(struct bp_array* bp_array, int idx, struct bp* bp);
struct bp*  get_bp(struct bp_array* bp_array, int idx);
char get_bp_chain_id(struct bp* bp, int side);
int get_bp_pos(struct bp* bp, int side);  
char* get_rna_sequence(struct rna_data* RNA);
char* set_chain_sequence(struct chain* RNA_chain);
char* get_chain_sequence(struct chain* RNA_chain);

#endif
