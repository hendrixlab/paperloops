import sys, os
from brewer2mpl import diverging, sequential
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def length(label):
    return sum([int(L) for L in label.split('-')])

def number(label):
    return len(label.split('-'))

def rotate(label):
    terms = label.split('-')
    new_terms = []
    for i in range(1,len(terms)):
        # effectively shifting over by 1 e.g. new_terms[i-1] = terms[i]
        new_terms.append(terms[i])
    new_terms.append(terms[0])
    rotation = '-'.join(new_terms)
    return rotation

def get_rep_label(label):
    if '-' in label:
        terms = label.split('-')
        rep_label = label
        for i in range(len(terms)):
            label = rotate(label)
            if label < rep_label:
                rep_label = label
    else:
        return label                
    return rep_label
    
def read_all_results_file(results_file):
    kB = 1.987204259e-3
    T = 310.15
    bmap = diverging.RdYlBu[11]
    colors = bmap.hex_colors
    F = open(results_file)
    data = {}
    Nloops = []
    labels = {}
    for line in F:
        if not line.startswith('motif'):
            line = line.replace(',','\t')
            if len(line.split()) == 3:
                label,C_eff,C_eff_err = line.split('\t')
                C_eff = float(C_eff)
                if length(label) not in colors:
                    color = colors[length(label) % 11]
                else:
                    print("could not find in colors:",length(label))
                    exit()
                n = number(label)
                if n == 1:
                    color = colors[7]
                if n not in data:
                    Nloops.append(n)
                    data[n] = {}
                    labels[n] = []
                #print("loading all data for n =",n);
                # append to lexicographically lowest rotation
                if 0.0 < C_eff and C_eff <= 1.0:
                    label = get_rep_label(label)
                    if label not in data[n]:
                        # label should be changed to rep_label if n>2
                        labels[n].append(label)
                        data[n][label] = []
                    deltaG = -kB * T * np.log(C_eff);            
                    data[n][label].append((deltaG,float(C_eff_err),2*float(C_eff)/float(C_eff_err),color))
    for n in labels:
        labels[n].sort(key=lambda x:length(x))
    return data,labels,Nloops


def read_results_file(results_file):
    F = open(results_file)
    kB = 1.987204259e-3
    T = 310.15
    plt_data = {}
    for line in F:
        if not line.startswith('motif'):
            line = line.replace(',','\t')
            label,C_eff,C_eff_err = line.split('\t')
            C_eff = float(C_eff)
            n = number(label)
            if n not in plt_data:
                plt_data[n] = []
            #plt_data[n].append((label,float(C_eff),float(C_eff_err)))
            if 0.0 < C_eff and C_eff <= 1.0:                
                deltaG = -kB * T * np.log(C_eff);            
                plt_data[n].append((label,deltaG,float(C_eff_err)))                
    return plt_data

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <C_eff_all_results.txt> <C_eff_results.txt>"

if len(sys.argv) != 3:
    print(usage)
    exit(0)

all_results_file = sys.argv[1]
results_file = sys.argv[2]

data,labels,Nloops = read_all_results_file(all_results_file)
plt_data = read_results_file(results_file)

fig_sizes = {1:(5,2), 2:(10,5), 3:(15,5), 4:(20,5), 5:(20,5), 6:(25,5), 7:(30,5), 8:(30,5)}

first = True
for n in Nloops:
    if n in fig_sizes:
        plt.figure(figsize=fig_sizes[n])
        # results file, just the main averages 
        x_labels,y_values,y_err = zip(*plt_data[n])    
        if n==1:
            # plot errorbar just for hairpin loops
            plt.plot(range(len(x_labels)),y_values,c='k',zorder=1,linewidth=1,label='inverse-error weighted average')
        for i in range(len(labels[n])):
            total = len(data[n][labels[n][i]])
            if(total):
                C_eff,C_eff_err,sizes,colors = zip(*data[n][labels[n][i]])
                if first:
                    plt.scatter(total*[i],C_eff,s=sizes,c=colors,alpha=1.0,edgecolors='black',linewidth=0.1,zorder=2,label='average for steptype rotation')
                    first = False
                else:
                    plt.scatter(total*[i],C_eff,s=sizes,c=colors,alpha=1.0,edgecolors='black',linewidth=0.1,zorder=2)   
            else:
                print(f"Found no examples in all_results: n={n} ",labels[n][i])
                exit()
        #plt.yscale('log')
        x_values = range(len(labels[n]))
        plt.xticks(x_values[::3],labels[n][::3],rotation=90, fontsize=8)
        if n == 1:
            plt.xlabel("loop length (nt)")
        else:
            plt.xlabel("loop motif length (nt)") 
        plt.ylabel(r'$\Delta G \; (kcal/mol)$')
        plt.legend(markerscale=4)
        plt.savefig('NN_Circularization_C_eff_n'+str(n)+'.pdf',bbox_inches='tight')
        plt.clf()
