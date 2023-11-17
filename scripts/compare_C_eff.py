import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os


def read_results_file(results_file):
    F = open(results_file)
    data = {}
    for line in F:
        if not line.startswith('motif'):
            line = line.replace(',','\t')
            label,C_eff,C_eff_err = line.split('\t')
            data[label] = float(C_eff)
    return data

usage = "Usage: " + sys.argv[0] + " <C_eff_results.txt> <comparison>"

if len(sys.argv) != 3:
    print(usage)
    exit(0)

results_file1 = sys.argv[1]
results_file2 = sys.argv[2]
    
results1 = read_results_file(results_file1)
results2 = read_results_file(results_file2)

x_values = []
y_values = []

labels = []
for label in results1:
    if label in results2:
        x_values.append(results1[label])
        y_values.append(results2[label])
        if results2[label] > 0.7:
            labels.append((label,(results1[label],results2[label])))

plt.scatter(x_values,y_values,s=3)
for item in labels:
    name,xy=item
    plt.annotate(name,xy)
plt.xlabel(results_file1)
plt.ylabel(results_file2)
#plt.xscale(log)
plt.yscale('symlog')
#plt.xlim(0,0.5)
#plt.ylim(0,10)
plt.savefig("C_eff_comparison_scatter.pdf")
