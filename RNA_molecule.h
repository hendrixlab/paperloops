/*
 * This file contains the definition of an interface for a linked list structure
 * of nucleotide nodes, each with a linked list of atom nodes. .
 */

#ifndef _RNA_MOLECULE
#define _RNA_MOLECULE

/*
 * Structure used to represent RNA molecule linked list.
 */
struct rna;

/*
 * Structure used to represent the nucleotide node.
 */
struct nuc;

/*
 * Structure used to represent the atom node.
 */
struct atom;


struct rna* create_rna();
struct nuc* create_nuc(int nuc_num, char* nuc_name);
struct atom* create_atom(int atom_num, char* atom_name, double* r);
void free_rna(struct rna* RNA);
void free_nuc(struct nuc* n);
void free_atom(struct atom* a);
int rna_is_empty(struct rna* RNA);
int nuc_is_empty(struct nuc* n);
void add_nuc(struct rna* RNA, struct nuc* n);
void add_atom(struct rna* RNA, struct atom* a);
int last_nuc_has_num(struct rna* RNA, int nuc_num);
struct nuc* get_nuc_by_num(struct rna* RNA, int nuc_num);
struct atom* get_atom_by_num(struct nuc* n, int atom_num);
int get_nuc_num(struct nuc* n);
char* get_nuc_name(struct nuc* n);
int get_atom_num(struct atom* a);
char* get_atom_name(struct atom* a);

#endif
