import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import sys, os, random
from scipy.optimize import curve_fit


def main():
    usage = "Usage: " + sys.argv[0] + " <C_eff_results.txt>"
    
    if len(sys.argv) != 2:
        print(usage)
        exit(0)
        
    results_file_list = sys.argv[1]
    
    combined = read_results_file_list(results_file_list)
    
    colors = ['#9e0142','#d53e4f','#f46d43','#fdae61','#abdda4','#66c2a5','#3288bd','#5e4fa2']

    data = {}
    max_n = 0
    for n in combined:
        if n > 1:
            data[n] = []
            for L in combined[n]:
                mean_dG = np.mean(combined[n][L])
                data[n].append((n,L,mean_dG,colors[n]))
                if n > max_n:
                    max_n = n
    n_order = [2,4,6,3,5,7]
    for n in n_order:
        # x_values are the lengths, y_values are the delta G values
        n_values,L_values,dG_values,cols = zip(*data[n])
        plt.scatter(L_values,dG_values,s=5,label='multiloops, n='+str(n),c=cols)
        # fitting a curve
        x = np.array(L_values)
        y = np.array(dG_values)
        n_values = np.array(n_values)
        L_values = np.array(L_values)
        X = (n_values,L_values)        
        p0 = (random.uniform(-5,0),random.uniform(0,5),random.uniform(0,5),1.5+noise(0.05),random.uniform(-5,0))
        if n == 2:
            p0 = (random.uniform(0,5),random.uniform(0,5),random.uniform(0,5),1.5+noise(0.05),random.uniform(0,5))
        if n > 3:
            p0 = (A_prev+noise(0.1*A_prev),B_prev+noise(0.1*B_prev),C_prev+noise(0.1*C_prev),D_prev+noise(0.1*D_prev),E_prev+noise(0.1*E_prev))
        print("initital:",n,p0)
        popt,pcov = curve_fit(f,X,y,p0,maxfev=10000,bounds=([-100.0,0.0,0.0,0.0,-100],100))
        A,B,C,D,E = popt
        eq_label = " %.2g*ln( %.2g*(L+n)^%.2g + %.2g*(n)^%.2g) + %.2g" % (A,B,D,C,D,E)
        plt.plot(x, f(X,*popt), color=cols[0], label=eq_label)
        print(n, eq_label)
        A_prev,B_prev,C_prev,D_prev,E_prev = A,B,C,D,E
    plt.xlabel('L, combined loop length (nt)')
    plt.ylabel(r'$\Delta$G loop circularization (kcal/mol)');
    plt.legend()
    plt.savefig("multiloop_equation_aalberts.pdf")

def noise(val):
    return random.uniform(-val,val)
    
def f(X,A,B,C,D,E):
    n,L = X
    kB = 1.987204259e-3
    T = 310.15
    return A * kB * T * np.log( B*np.power(L+n,D) + C*np.power(n,D) ) + E
    
def length(label):
    return sum([int(L) for L in label.split('-')])

def count(label):
    return len(label.split('-'))

#########
# INPUT #
#########

def read_results_file_list(results_file_list):
    F = open(results_file_list)
    combined = {}
    for line in F:
        results_file = line.strip()
        results = read_results_file(results_file)
        for n in results:
            if n not in combined:
                combined[n] = {}
            for L in results[n]:
                if L not in combined[n]:
                    combined[n][L] = []
                combined[n][L].extend(results[n][L])                
    return combined

def read_results_file(results_file):
    F = open(results_file)
    kB = 1.987204259e-3
    T = 310.15
    results = {}
    for line in F:
        if not line.startswith('motif'):
            line = line.replace(',','\t')
            label,C_eff,C_eff_err = line.split('\t')
            C_eff = float(C_eff)
            if C_eff > 0.0:
                deltaG = -kB * T * np.log(C_eff);
                n = count(label)
                if n not in results:
                    results[n] = {}
                L = length(label)
                if L not in results[n]:
                    results[n][L] = []
                results[n][L].append(deltaG)
    return results


########
# MAIN #
########

main()
