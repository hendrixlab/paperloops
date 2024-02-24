SGE_Batch -c "./bin/paperloops -M loops_n7.txt -P data/4ybb_23S.pdb -S data/4ybb_23S.pdb.stems.txt -x large3_rep1" -r large3_rep1 -q cascade
sleep 20
SGE_Batch -c "./bin/paperloops -M loops_n7.txt -P data/4ybb_23S.pdb -S data/4ybb_23S.pdb.stems.txt -x large3_rep2" -r large3_rep2 -q cascade
sleep 20
SGE_Batch -c "./bin/paperloops -M loops_n7.txt -P data/4ybb_23S.pdb -S data/4ybb_23S.pdb.stems.txt -x large3_rep3" -r large3_rep3 -q cascade
sleep 20
SGE_Batch -c "./bin/paperloops -M loops_n7.txt -P data/4ybb_23S.pdb -S data/4ybb_23S.pdb.stems.txt -x large3_rep4" -r large3_rep4 -q cascade
sleep 20
SGE_Batch -c "./bin/paperloops -M loops_n7.txt -P data/4ybb_23S.pdb -S data/4ybb_23S.pdb.stems.txt -x large3_rep5" -r large3_rep5 -q cascade
