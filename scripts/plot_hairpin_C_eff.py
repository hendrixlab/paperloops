import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os

def read_results_file(results_file):
    F = open(results_file)
    data = []
    for line in F:
        if not line.startswith('motif'):
            line = line.replace(',','\t')
            label,C_eff,C_eff_err = line.split('\t')
            if '-' not in label:
                data.append((label,float(C_eff),float(C_eff_err)))
    return data

usage = "Usage: " + sys.argv[0] + " <C_eff_results.txt>"

if len(sys.argv) != 2:
    print(usage)
    exit(0)

results_file = sys.argv[1]

plot_data = read_results_file(results_file)
plt_data.sort(key=lambda x:x[0])
x_values,y_values,y_err = zip(*plt_data)

plt.errorbar(x_values,y_values,y_err)
plt.xlabel("Hairpin loop length (nt)")
plt.ylabel("Estimated effective concentration (M)")
plt.savefig('hairpin_loops_C_eff.pdf')
