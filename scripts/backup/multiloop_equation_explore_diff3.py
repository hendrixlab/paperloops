import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from scipy.optimize import curve_fit
from matplotlib.colors import LinearSegmentedColormap
import random


def main():
    usage = "Usage: " + sys.argv[0] + " <C_eff_results.txt>"
    
    if len(sys.argv) != 2:
        print(usage)
        exit(0)
        
    results_file_list = sys.argv[1]
    
    combined = read_results_file_list(results_file_list)
    
    colors = ['#9e0142','#d53e4f','#f46d43','#fdae61','#abdda4','#66c2a5','#3288bd','#5e4fa2']

    kB = 1.987204259e-3
    T = 310.15
    data_short = {}
    avg_all = {}
    data = {}
    max_n = 0
    for n in combined:
        if n > 1:
            data[n] = []
            for L in combined[n]:
                mean_dG = np.mean(combined[n][L])                
                if L > 0:
                    if L not in avg_all:
                        avg_all[L] = []
                    avg_all[L].append(mean_dG)
                    data[n].append((n,L,mean_dG,colors[n]))

    # averages for all loops
    L_values = []
    dG_values = []
    for L in avg_all:
        avg_dG = np.mean(avg_all[L])
        L_values.append(L)
        dG_values.append(avg_dG)
    L_values = np.array(L_values)
    dG_values = np.array(dG_values)
    #p0 = (A,B,D)
    popt_ideal,pcov = curve_fit(f_ideal,L_values,dG_values,maxfev=100000,bounds=([-100,0.0,-100.0],100))
    (A,B,C) = popt_ideal
    eq_label = r'$\Delta G_{avg}$ ='+str(round(A,2))+'*kB*T*log(x + ' + str(round(B,2)) + ') + '+ str(round(C,2))
    plt.scatter(L_values,dG_values,s=5,label='Average',c='gray')
    plt.plot(L_values, f_ideal(L_values,*popt_ideal), color='gray', label=eq_label)

    # plot the full data set 
    params = []
    popt_dict = {}    
    keys = sorted(list(data.keys()))
    for n in keys:
        # x_values are the lengths, y_values are the delta G values
        n_values,L_values,dG_values,cols = zip(*data[n])
        # plotting # fitting
        L_values = np.array(L_values)
        y_values_diff = np.array(dG_values) - f_ideal(L_values,*popt_ideal)
        y_values = np.array(dG_values)
        plt.scatter(L_values,y_values,s=5,label='multiloops, n='+str(n),c=cols)
        '''
        p0 = (3.1411,3.4583737,5.42421)
        if n % 2 == 0:
            p0 = (-3.1411,3.4583737,5.42421)
        '''
        p0 = (3.1411,3.4583737)
        if n % 2 == 0:
            p0 = (-3.1411,3.4583737)
        popt,pcov = curve_fit(g,L_values,y_values_diff,p0,maxfev=100000,bounds=([-100,0.0],100))
        A,B,C = popt_ideal
        D,E = popt
        popt_combined = (A,B,C,D,E)
        eq_label = r'$\Delta\Delta$G=' + str(round(A,2)) + '*kB*T*log(x + ' + str(round(B,2)) + ') + ' + str(round(D,2))+'/(x^2 +  ' + str(round(E,2)) + ') '
        #eq_label = r'$\Delta$G='+
        plt.plot(L_values, f_ideal_g(L_values,*popt_combined), color=cols[0], label=eq_label)
    plt.xlabel('L, combined loop length (nt)')
    plt.ylabel(r'$\Delta$G loop circularization (kcal/mol)');
    plt.legend(fontsize=8)
    plt.savefig("multiloop_equation_explore_diff3.pdf")
    plt.clf()

def f_ideal(x,A,B,C):
    kB = 1.987204259e-3
    T = 310.15
    return A*kB*T*np.log(x + B) + C

def g(x,D,E):
    return D/(x**2 + E)

def f_ideal_g(x,A,B,C,D,E):
    kB = 1.987204259e-3
    T = 310.15
    return A*kB*T*np.log(x+B) + C + D/(x**2 + E)

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
