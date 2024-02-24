import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from scipy.optimize import curve_fit
from matplotlib.colors import LinearSegmentedColormap


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
    data = {}
    max_n = 0
    const = {}
    for n in combined:
        if n > 1:
            data[n] = []
            for L in combined[n]:
                mean_dG = np.mean(combined[n][L])
                if L == 0:
                    const[n] = mean_dG/(kB*T) # this is the y-intercept
                    print(n,const[n]);
                data[n].append((n,L,mean_dG,colors[n]))
                if n > max_n:
                    max_n = n

    for n in data:
        # x_values are the lengths, y_values are the delta G values
        n_values,L_values,dG_values,cols = zip(*data[n])
        plt.scatter(L_values,dG_values,s=5,label='multiloops, n='+str(n),c=cols)
        # fitting a curve
        x = np.array(L_values)
        y = np.array(dG_values)
        p0 = (const[n],1.0137,3.03070,1.112149,0.10101)
        popt,pcov = curve_fit(f,x,y,p0)
        A,B,C,D,E = popt
        eq_label = r'$\Delta$G='+str(round(A,2))+'*log('+str(round(B,2))+' + ' + str(round(C,2))+'*x^'+str(round(D,2)) + ') + ' + str(round(E,2))
        #if B < 0:
        #    eq_label = r'$\Delta$G='+str(round(A,2))+' - '+str(round(abs(B),2))+'*x'            
        plt.plot(x, f(x,*popt), color=cols[0], label=eq_label)
        print(n,*popt)


    plt.xlabel('L, combined loop length (nt)')
    plt.ylabel(r'$\Delta$G loop circularization (kcal/mol)');
    plt.legend()
    plt.savefig("multiloop_equation_logarithmic.pdf")

def f(x,A,B,C,D,E):
    kB = 1.987204259e-3
    T = 310.15
    return A*kB*T*np.log(B + C*x**D) + E
    
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
