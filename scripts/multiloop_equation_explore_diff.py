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
    data_long = {}
    data = {}
    avg_long = {}
    avg_short = {}
    max_n = 0
    for n in combined:
        if n > 1:
            data_short[n] = []
            data_long[n] = []
            data[n] = []
            for L in combined[n]:
                mean_dG = np.mean(combined[n][L])                
                if L > 8:
                    # the average as a function of L for large values
                    if L not in avg_long:
                        avg_long[L] = []
                    avg_long[L].append(mean_dG)
                    data_long[n].append((n,L,mean_dG,colors[n]))
                elif L > 0: # should be from 1 and 8 inclusive
                    if L not in avg_short:
                        avg_short[L] = []
                    avg_short[L].append(mean_dG)
                    data_short[n].append((n,L,mean_dG,colors[n]))
                if n > max_n:
                    max_n = n
                if L > 0:
                    data[n].append((n,L,mean_dG,colors[n]))

    # compute some averages for short loops:
    L_values = []
    dG_values = []
    for L in avg_short:
        avg_dG = np.mean(avg_short[L])
        L_values.append(L)
        dG_values.append(avg_dG)
    L_values = np.array(L_values)
    dG_values = np.array(dG_values)
    #p0 = (A,B,D)
    popt,pcov = curve_fit(f2,L_values,dG_values,maxfev=100000,bounds=([-100,0.0,-100.0],100))
    (A,B,D) = popt
    eq_label = r'$\Delta G_{avg,short}$ ='+str(round(A,2))+'*kB*T*log(x + ' + str(round(B,2)) + ') + '+ str(round(D,2))
    plt.scatter(L_values,dG_values,s=5,label='Average',c='gray')
    plt.plot(L_values, f2(L_values,*popt), color='gray', label=eq_label)

    # averages for the long loops
    L_values = []
    dG_values = []
    for L in avg_long:
        avg_dG = np.mean(avg_long[L])
        L_values.append(L)
        dG_values.append(avg_dG)
    L_values = np.array(L_values)
    dG_values = np.array(dG_values)
    #p0 = (A,B,D)
    popt_long,pcov = curve_fit(f2,L_values,dG_values,maxfev=100000,bounds=([-100,0.0,-100.0],100))
    (A,B,D) = popt_long
    eq_label = r'$\Delta G_{avg,long}$ ='+str(round(A,2))+'*kB*T*log(x + ' + str(round(B,2)) + ') + '+ str(round(D,2))
    plt.scatter(L_values,dG_values,s=5,label='Average',c='gray')
    plt.plot(L_values, f2(L_values,*popt_long), color='gray', label=eq_label)
    plt.xlabel('L, combined loop length (nt)')
    plt.ylabel(r'$\Delta$G loop circularization (kcal/mol)');
    plt.legend(fontsize=8)
    plt.savefig("multiloop_equation_confirm_diff.pdf")
    plt.clf()

    # plot the full data set 
    params = []
    popt_dict = {}    
    keys = sorted(list(data_short.keys()))
    for n in keys:
        # x_values are the lengths, y_values are the delta G values
        n_values,L_values,dG_values,cols = zip(*data_short[n])
        # plotting
        if data_long[n]:
            n_values_long,L_values_long,dG_values_long,cols_long = zip(*data_long[n])
            L_values = L_values + L_values_long
            dG_values = dG_values + dG_values_long
            cols = cols + cols_long
        # fitting
        L_values = np.array(L_values)
        y_values = np.array(dG_values) - f2(L_values,*popt_long)
        plt.scatter(L_values,y_values,s=5,label='multiloops, n='+str(n),c=cols)
        '''
        p0 = (3.1411,3.4583737,5.42421)
        if n % 2 == 0:
            p0 = (-3.1411,3.4583737,5.42421)
        '''
        p0 = (3.1411,3.4583737)
        if n % 2 == 0:
            p0 = (-3.1411,3.4583737)
        popt,pcov = curve_fit(g,L_values,y_values,p0,bounds=([-100,0.0],100))
        A,B = popt
        #A,B,C = popt        
        #popt_dict[n] = (A,B,C,D)
        #params.append((n,A,B,C,D))
        #eq_label = r'$\Delta\Delta$G='+str(round(A,2))+'/(x +  ' + str(round(B,2)) + ') - ' + str(round(C,2))
        eq_label = r'$\Delta\Delta$G='+str(round(A,2))+'/(x +  ' + str(round(B,2)) + ') '
        #eq_label = r'$\Delta$G='+str(round(A,2))+'*kB*T*log(x + ' + str(round(B,2)) + ') - ' + str(round(C,2)) + '*x + ' + str(round(D,2))
        plt.plot(L_values, g(L_values,*popt), color=cols[0], label=eq_label)
    plt.xlabel('L, combined loop length (nt)')
    plt.ylabel(r'$\Delta$G loop circularization (kcal/mol)');
    plt.legend(fontsize=8)
    plt.savefig("multiloop_equation_explore_diff.pdf")
    plt.clf()

def f(x,A,B,C,D):
    #n,L = x
    kB = 1.987204259e-3
    T = 310.15
    #print("f: ",x,A,B,C,D)
    return A*kB*T*np.log(x + B) - C*x  + D

def f2(x,A,B,D):
    #n,L = x
    kB = 1.987204259e-3
    T = 310.15
    #print("f: ",x,A,B,C,D)
    return A*kB*T*np.log(x + B) + D

def f3(x,A,B,C,D):
    #n,L = x
    kB = 1.987204259e-3
    T = 310.15
    #print("f: ",x,A,B,C,D)
    return A*kB*T*np.log(x + B) - C*kB*T*np.log(x)+ D

'''

def f(x,A,B,C):
    #kB = 1.987204259e-3
    #T = 310.15
    return A - B*np.log(x+C)

def g(x,A,B,C):
    return A/(x**2 + B) + C

'''
def g(x,A,B):
    return A/(x**2 + B)

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
