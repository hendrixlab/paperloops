CC=cc -std=c99 -g

paperloops: paperloops.c rna_data.o simulation.o transform.o rna_io.o motifs.o vectorfun.o 
	$(CC) -fopenmp paperloops.c rna_data.o rna_io.o simulation.o transform.o motifs.o vectorfun.o -o ../bin/paperloops -lm

simulation.o: simulation.c simulation.h transform.c transform.h transform.o vectorfun.c vectorfun.h vectorfun.o
	$(CC) -c simulation.c vectorfun.o transform.o -lm

transform.o: transform.c transform.h rna_data.o motifs.o vectorfun.o
	$(CC) -c transform.c rna_data.o vectorfun.o motifs.o -o transform.o -lm

rna_io.o: rna_io.c rna_io.h rna_data.o vectorfun.o
	$(CC) -c rna_io.c rna_data.o transform.o vectorfun.o -o rna_io.o -lm

rna_data.o: rna_data.c rna_data.h vectorfun.o motifs.o
	$(CC) -c rna_data.c vectorfun.o motifs.o -o rna_data.o -lm

motifs.o: motifs.c motifs.h
	$(CC) -c motifs.c -o motifs.o

vectorfun.o: vectorfun.c vectorfun.h
	$(CC) -c vectorfun.c -o vectorfun.o

clean:
	rm -f simulation.o transform.o rna_io.o rna_data.o vectorfun.o motifs.o
	rm -f *~
