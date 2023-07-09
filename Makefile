CC=cc --std=c99 -g

read_RNA_PDB: read_RNA_PDB.c RNA_molecule.o
	$(CC) read_RNA_PDB.c RNA_molecule.o -o read_RNA_PDB

RNA_molecule.o: RNA_molecule.c RNA_molecule.h
	$(CC) -c RNA_molecule.c

clean:
	rm -f read_RNA_PDB RNA_molecule.o 
	rm -rf read_RNA_PDB.dSYM
