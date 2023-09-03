CC=cc --std=c99 -g

test: test.c rna_data.o transform.o rna_pdb.o vectorfun.o
	$(CC) test.c rna_data.o rna_pdb.o transform.o vectorfun.o -o test -lm

transform.o: transform.c transform.h rna_data.o vectorfun.o
	$(CC) -c transform.c rna_data.o vectorfun.o -o transform.o -lm

rna_pdb.o: rna_pdb.c rna_pdb.h rna_data.o vectorfun.o
	$(CC) -c rna_pdb.c rna_data.o transform.o vectorfun.o -o rna_pdb.o -lm

rna_data.o: rna_data.c rna_data.h vectorfun.o
	$(CC) -c rna_data.c vectorfun.o -o rna_data.o -lm

vectorfun.o: vectorfun.c vectorfun.h
	$(CC) -c vectorfun.c -o vectorfun.o

clean:
	rm -f rna_pdb.o rna_data.o transform.o vectorfun.o
	rm -rf read_rna_pdb.dSYM
	rm -f test
	rm -f *sugar.txt *backbone.txt
	rm -f *~
