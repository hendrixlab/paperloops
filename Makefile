CC=cc --std=c99 -g

test: test.c RNA_molecule.o RNA_PDB.o
	$(CC) test.c RNA_molecule.o RNA_PDB.o -o test

RNA_PDB.o: RNA_PDB.c RNA_PDB.h
	$(CC) -c RNA_PDB.c 

RNA_molecule.o: RNA_molecule.c RNA_molecule.h
	$(CC) -c RNA_molecule.c 

clean:
	rm -f RNA_PDB.o RNA_molecule.o 
	rm -rf test.dSYM
	rm -f *~
