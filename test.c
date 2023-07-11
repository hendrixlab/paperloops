#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "RNA_molecule.h"
#include "RNA_PDB.h"

int main(int argc, char *argv[]) {
  if(argc == 2) {
    printf("input file name: %s\n",argv[1]);
  } else {
    printf("Usage: %s <PDB file name>\n",argv[0]);
    exit(1);
  }
  char* file_name = argv[1];
  int model_count;
  struct rna** models = read_pdb_file(file_name,&model_count);

  struct rna* RNA = models[0];
  struct nuc* n = get_first_nuc(RNA);
  int FIRST = 1;
  while(n) {
    printf("%d - %s:",get_nuc_num(n),get_nuc_name(n));
    /*
    if(!FIRST) {
      struct atom* prev_a = get_prev_backbone_atom(n);
    }
    */
    struct atom* a = get_first_backbone_atom(n);
    while(a) {
      printf(" - %s",get_atom_name(a));
      a = next_atom(a);      
    }    
    printf("\n");
    n = next_nuc(n);
  }
  char* sequence = get_rna_sequence(RNA);
  printf("%s\n",sequence);

  /*
  // this is some example code to confirm that the data is stored properly
  struct rna* RNA = models[0];
  struct nuc* n = get_nuc_by_num(RNA,1);
  struct atom* a = get_atom_by_num(n,7);
  printf("nuc num: %d nuc name: %s\n",get_nuc_num(n),get_nuc_name(n));
  printf("atom num: %d, atom name: %s, atom type: %s\n",get_atom_num(a),get_atom_name(a),get_atom_type(a));

  // print the RNA sequence:
  char* sequence = get_rna_sequence(RNA);
  printf("%s\n",sequence);
  */
}
