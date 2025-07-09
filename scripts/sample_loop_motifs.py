import sys, os, re
import random

def main():
    usage = "Usage: " + sys.argv[0] + " <loop motifs> <min_n> <min_L> <max_L>"
    if len(sys.argv) != 6:
        print(usage)
        exit()

    loop_motifs_file = sys.argv[1]
    min_n = int(sys.argv[2])
    min_L = int(sys.argv[3])
    max_L = int(sys.argv[4])
    size = int(sys.argv[5])

    loop_motifs = read_loop_motifs(loop_motifs_file)

    for n in loop_motifs:
        for L in loop_motifs[n]:
            if min_L <= L and L <= max_L and min_n <= n:
                #print(f"# n={n} L={L}")
                if len(loop_motifs[n][L]) > size:
                    # sample motifs:                
                    motif_sample = random.sample(loop_motifs[n][L],k=size)
                    for loop_motif in motif_sample:
                        print(loop_motif)
                else:
                    for loop_motif in loop_motifs[n][L]:
                        print(loop_motif)


#############
# FUNCTIONS #
#############

def sample_motifs(motif_list,size):
    return random.sample(motif_list,k=size)

def read_loop_motifs(loop_motifs_file):
    F = open(loop_motifs_file)
    loop_motifs = {}
    n = 0
    L = 0
    for line in F:
        if line.startswith('#'):
            # Npred(n=1,L=0)=1, Nobs=1
            match = re.search('Npred\(n=(\d+),L=(\d+)\)',line)
            n = int(match.group(1))
            L = int(match.group(2))
            if n not in loop_motifs:
                loop_motifs[n] = {}
            if L not in loop_motifs[n]:
                loop_motifs[n][L] = []
        else:
            # variable number of collumns
            motif_info = line.strip().split()
            number = motif_info[0]
            lword = motif_info[1] # currenlty storing only the lowest motif
            loop_motifs[n][L].append(lword)
    return loop_motifs

main()
