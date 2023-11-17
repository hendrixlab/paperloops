import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from brewer2mpl import diverging

def length(label):
    return sum([int(L) for L in label.split('-')])

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
    terms = label.split('-')
    rep_label = label
    for i in range(len(terms)):
        label = rotate(label)
        if label < rep_label:
            rep_label = label
    return rep_label
    
def read_results_file(results_file):
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
                if length(label) not in colors:
                    color = colors[length(label)]
                else:
                    print("could not find in colors:",length(label))
                    exit()
                n = len(label.split('-'))
                if n not in data:
                    Nloops.append(n)
                    data[n] = {}
                    labels[n] = []
                if label != '0':                    
                    if n > 1:
                        # append to lexicographically lowest rotation
                        label = get_rep_label(label)
                    if label not in data[n]:
                        # label should be changed to rep_label if n>2
                        labels[n].append(label)
                        data[n][label] = []
                    else:
                        if float(C_eff_err) != 0.0:
                            data[n][label].append((float(C_eff),float(C_eff_err),float(C_eff)/float(C_eff_err),color))
                        else:
                            data[n][label].append((float(C_eff),float(C_eff_err),0.0,color))
    return data,labels,Nloops

usage = "Usage: " + sys.argv[0] + " <C_eff_all_results.txt>"

if len(sys.argv) != 2:
    print(usage)
    exit(0)

results_file = sys.argv[1]

data,labels,Nloops = read_results_file(results_file)

fig,axs = plt.subplots(nrows=len(Nloops), figsize=(20,15))

for n in Nloops:
    for i in range(len(labels[n])):
        total = len(data[n][labels[n][i]])
        C_eff,C_eff_err,sizes,colors = zip(*data[n][labels[n][i]])
        axs[n-1].scatter(total*[i],C_eff,s=sizes,c=colors,alpha=0.5,edgecolors=None)
    axs[n-1].set_yscale('log')
    axs[n-1].set_xticks(range(len(labels[n])))
    axs[n-1].set_xticklabels(labels[n],rotation=90, fontsize=8)
    axs[n-1].set_xlabel("loop types")

fig.supylabel("Estimated effective concentration (M)")

plt.savefig('Nearest_Neighbor_Circularization_C_eff.pdf',bbox_inches='tight')
