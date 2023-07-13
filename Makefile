CC=cc --std=c99 -g

test: test.c RNA_molecule.o RNA_PDB.o
	$(CC) test.c RNA_molecule.o RNA_PDB.o vectorfun.o -o test -lm

RNA_PDB.o: RNA_PDB.c RNA_PDB.h RNA_molecule.o
	$(CC) -c RNA_PDB.c RNA_molecule.o -o RNA_PDB.o -lm

RNA_molecule.o: RNA_molecule.c RNA_molecule.h vectorfun.o
	$(CC) -c RNA_molecule.c vectorfun.o -o RNA_molecule.o -lm

vectorfun.o: vectorfun.c vectorfun.h
	$(CC) -c vectorfun.c -o vectorfun.o

clean:
	rm -f RNA_PDB.o RNA_molecule.o vectorfun.o
	rm -rf read_RNA_PDB.dSYM
	rm -f test
	rm -f *~
