CC=cc --std=c99 -g
all: test1_build_libraries test2_list_motifs test3_motif_to_steptypes test4_sample_C_eff test5_sample_C_eff

test5_sample_C_eff: test5_sample_C_eff.c rna_data.o rna_io.o transform.o motifs.o simulation.o vectorfun.o
	$(CC) test5_sample_C_eff.c rna_data.o rna_io.o simulation.o transform.o motifs.o vectorfun.o -o test5_sample_C_eff -lm

test4_sample_C_eff: test4_sample_C_eff.c rna_data.o rna_io.o transform.o motifs.o simulation.o vectorfun.o
	$(CC) test4_sample_C_eff.c rna_data.o rna_io.o simulation.o transform.o motifs.o vectorfun.o -o test4_sample_C_eff -lm

test3_motif_to_steptypes: test3_motif_to_steptypes.c rna_data.o rna_io.o transform.o motifs.o simulation.o vectorfun.o
	$(CC) test3_motif_to_steptypes.c rna_data.o rna_io.o simulation.o transform.o motifs.o vectorfun.o -o test3_motif_to_steptypes -lm

test2_list_motifs: test2_list_motifs.c rna_data.o transform.o rna_io.o motifs.o vectorfun.o
	$(CC) test2_list_motifs.c rna_data.o rna_io.o transform.o motifs.o vectorfun.o -o test2_list_motifs -lm

test1_build_libraries: test1_build_libraries.c rna_data.o transform.o rna_io.o motifs.o vectorfun.o 
	$(CC) test1_build_libraries.c rna_data.o rna_io.o transform.o motifs.o vectorfun.o -o test1_build_libraries -lm

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
	rm -f rna_io.o rna_data.o transform.o vectorfun.o motifs.o
	rm -f test1_build_libraries test2_list_motifs test3_motif_to_steptypes test4_sample_C_eff test5_sample_C_eff
	rm -f *sugar.txt *backbone.txt
	rm -f *~
