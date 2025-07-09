import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from scipy.optimize import curve_fit
from matplotlib.colors import LinearSegmentedColormap
import random
from sklearn.metrics import r2_score

def main():
    usage = "Usage: " + sys.argv[0] + " <collision file 1> <collision file 2>"
    if len(sys.argv) != 3:
        print(usage)
        exit()

    data_file1 = sys.argv[1]
    data_file2 = sys.argv[2]
    
    data1 = read_data(data_file1)
    data2 = read_data(data_file2)

    plt.hist(data1)
    plt.hist(data2)

    plt.savefig('test.pdf')
    
#############
# FUNCTIONS #
#############

def read_data(data_file):
    F = open(data_file)
    data = []
    for line in F:
        if not line.startswith('motif'):
            motif,steps,p,count,total,fwd_collide,rev_collide,N_samples = line.strip().split('\t')
            data.append(float(fwd_collide))
    return data

main()
