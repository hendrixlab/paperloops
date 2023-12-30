import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from brewer2mpl import diverging
from matplotlib.cm import ScalarMappable
from scipy.optimize import curve_fit


def main():
    usage = "Usage: " + sys.argv[0] + " <C_eff_results.txt> <comparison>"
    
    if len(sys.argv) != 3:
        print(usage)
        exit(0)
        
    results_file1 = sys.argv[1]
    results_file2 = sys.argv[2]
    
    results1,max_length = read_results_file(results_file1)
    results2,max_length2 = read_results_file(results_file2)
    
    bmap = diverging.RdYlBu[min(11,max_length+1)]
    colors = bmap.hex_colors
    
    data = []
    error = {}
    for label in results1:
        if label in results2:
            L = length(label)
            if not np.isnan(results2[label][0]):
                data.append((results1[label][0],results2[label][0],label,colors[L % 11]))
                error[label] = results1[label][1]
    data.sort(key=lambda x:x[0])
    x_values,y_values,labels,colors = zip(*data)
    plt.plot(list(range(10)),list(range(10)),c='k',alpha=0.3,label='y=x')
    plt.scatter(x_values,y_values,s=3,label='N_loops=5, L_max=5',color='b')
    
    x = np.array(x_values)
    y = np.array(y_values)
    popt,pcov = curve_fit(f,x,y)
    x = np.unique(x)
    plt.plot(x, f(x,*popt), label="y=%.3fx+%.3f" % tuple(popt), color='r')
    
    plt.xlabel(results_file1)
    plt.ylabel(results_file2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("paperloops")
    plt.ylabel("toyfold3")
    
    for item in data:
        x_value,y_value,label,color = item
        if abs(x_value - y_value)/x_value >= 0.25:
            xy = (x_value,y_value)
            plt.annotate(label,xy)


    plt.legend()
    #plt.xlim(0,0.5)
    #plt.ylim(0,10)
    #cm = bmap.get_mpl_colormap(N=max_length+1)
    #plt.colorbar(ScalarMappable(cmap=cm))
    plt.savefig("C_eff_comparison_scatter.pdf")
    
def length(label):
    return sum([int(L) for L in label.split('-')])

def read_results_file(results_file):
    F = open(results_file)
    data = {}
    max_length = 0
    for line in F:
        if not line.startswith('motif'):
            line = line.replace(',','\t')
            label,C_eff,C_eff_err = line.split('\t')
            data[label] = (float(C_eff),float(C_eff_err))
            L = length(label)
            if L > max_length:
                max_length = L
    return data,max_length

def f(x,m,b):
    return x*m + b

########
# MAIN #
########

main()
