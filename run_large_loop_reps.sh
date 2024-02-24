SGE_Batch -c "./bin/paperloops -M large_loops.txt -P data/4ybb_23S.pdb -S data/4ybb_23S.pdb.stems.txt -x large2_rep1" -r large2_rep1 -q cascade
sleep 20
SGE_Batch -c "./bin/paperloops -M large_loops.txt -P data/4ybb_23S.pdb -S data/4ybb_23S.pdb.stems.txt -x large2_rep2" -r large2_rep2 -q cascade
sleep 20
SGE_Batch -c "./bin/paperloops -M large_loops.txt -P data/4ybb_23S.pdb -S data/4ybb_23S.pdb.stems.txt -x large2_rep3" -r large2_rep3 -q cascade
sleep 20
SGE_Batch -c "./bin/paperloops -M large_loops.txt -P data/4ybb_23S.pdb -S data/4ybb_23S.pdb.stems.txt -x large2_rep4" -r large2_rep4 -q cascade
sleep 20
SGE_Batch -c "./bin/paperloops -M large_loops.txt -P data/4ybb_23S.pdb -S data/4ybb_23S.pdb.stems.txt -x large2_rep5" -r large2_rep5 -q cascade
