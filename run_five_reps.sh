rm -rf rep1
rm -rf rep2
rm -rf rep3
rm -rf rep4
rm -rf rep5
SGE_Batch -c "./test5_sample_C_eff data/4ybb_23S.pdb data/4ybb_23S.pdb.stems.txt rep1 > log.rep1.txt" -r rep1 -q cascade
sleep 10
SGE_Batch -c "./test5_sample_C_eff data/4ybb_23S.pdb data/4ybb_23S.pdb.stems.txt rep2 > log.rep2.txt" -r rep2 -q cascade
sleep 10
SGE_Batch -c "./test5_sample_C_eff data/4ybb_23S.pdb data/4ybb_23S.pdb.stems.txt rep3 > log.rep3.txt" -r rep3 -q cascade
sleep 10
SGE_Batch -c "./test5_sample_C_eff data/4ybb_23S.pdb data/4ybb_23S.pdb.stems.txt rep4 > log.rep4.txt" -r rep4 -q cascade
sleep 10
SGE_Batch -c "./test5_sample_C_eff data/4ybb_23S.pdb data/4ybb_23S.pdb.stems.txt rep5 > log.rep5.txt" -r rep5 -q cascade
