CC = cc
CFLAGS = -std=c99 -O3 -march=native -ffast-math -fopenmp

paperloops: paperloops.c rna_data.o simulation.o transform.o rna_io.o motifs.o vectorfun.o
	$(CC) $(CFLAGS) paperloops.c rna_data.o rna_io.o simulation.o transform.o motifs.o vectorfun.o -o ../bin/paperloops -lm

simulation.o: simulation.c simulation.h transform.h vectorfun.h
	$(CC) $(CFLAGS) -c simulation.c -o simulation.o -lm

transform.o: transform.c transform.h rna_data.h motifs.h vectorfun.h
	$(CC) $(CFLAGS) -c transform.c -o transform.o -lm

rna_io.o: rna_io.c rna_io.h rna_data.h vectorfun.h
	$(CC) $(CFLAGS) -c rna_io.c -o rna_io.o -lm

rna_data.o: rna_data.c rna_data.h vectorfun.h motifs.h
	$(CC) $(CFLAGS) -c rna_data.c -o rna_data.o -lm

motifs.o: motifs.c motifs.h
	$(CC) $(CFLAGS) -c motifs.c -o motifs.o

vectorfun.o: vectorfun.c vectorfun.h
	$(CC) $(CFLAGS) -c vectorfun.c -o vectorfun.o

clean:
	rm -f simulation.o transform.o rna_io.o rna_data.o vectorfun.o motifs.o
	rm -f *~
