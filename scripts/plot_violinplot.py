import sys, os, re
from brewer2mpl import diverging
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from math import log

########
# MAIN #
########

def main():
    usage = "Usage: " + sys.argv[0] + " <paperloops outpt> <toyfold3 output>"

    if len(sys.argv) != 3:
        print(usage)
        exit(0)

    paperloops_output_file = sys.argv[1]
    toyfold3_output_file = sys.argv[2]

    data1 = read_paperloops(paperloops_output_file)
    data2 = read_toyfold3(toyfold3_output_file)

    #paperloops = [data1['1'],data1['2'],data1['3'],data1['4'],data1['5'],data1['6'],data1['7'],data1['8'],data1['9'],data1['10']]
    #toyfold3 =  [data2['1'],data2['2'],data2['3'],data2['4'],data2['5'],data1['6'],data1['7'],data1['8'],data1['9'],data1['10']]
    paperloops = [data1['0-0'],data1['0-1'],data1['1-0']]
    toyfold3 =  [data2['0-0'],data2['0-1'],data2['1-0']]
    data = paperloops + toyfold3
    
    fig = plt.figure(figsize=(10,5))
    #plt.yscale('log')
    plt.violinplot(data)
    plt.xticks(list(range(1,6+1)),['0-0','0-1','1-0','0-0','0-1','1-0'])
    plt.figtext(0.3, 0.01, "paperloops", wrap=True, horizontalalignment='center', fontsize=12)
    plt.figtext(0.75, 0.01, "toyfold3", wrap=True, horizontalalignment='center', fontsize=12)
    #plt.ylabel('C_eff')
    #plt.savefig('violinplot_C_eff.pdf')
    #plt.ylabel('w')
    #plt.savefig('violinplot_w.pdf')
    #plt.ylabel('v')
    #plt.savefig('violinplot_v.pdf')
    plt.ylabel('-log(p)')
    plt.savefig('violinplot_p.pdf')
    #plt.ylabel('KDE bandwidth')
    #plt.savefig('violinplot_bw.pdf')
    #plt.ylabel('denominator of C_eff calculation')
    #plt.savefig('violinplot_den.pdf')
        
def read_paperloops(paperloops_output_file):
    data = {}
    with open(paperloops_output_file) as F:
        L = ""
        for line in F:
            if line.startswith('N_pred'):
                # N_pred(n=1,L=0)=1
                match = re.search('N_pred\(n=\d+,L=(\d+)\)',line)
                L = match.group(1)
                line = next(F)
                motif = line.strip()
                if motif not in data:
                    data[motif] = []
            elif re.search('\d\-\d',line):
                motif = line.strip()
                if motif not in data:
                    data[motif] = []
                print(line)
            else:
                print(f"storing data for paperloops: '{motif}'")
                x = line.strip()
                if float(x) != 0.0: data[motif].append(-log(float(x)))
                #if float(x) != 0.0: data[motif].append(float(x))
    return data
                
def read_toyfold3(toyfold3_output_file):
    data = {}
    with open(toyfold3_output_file,'r+') as F:
        motif = ""
        for line in F:
            if line.startswith('motif'):
                # motif=1
                match = re.search('motif=(.*)',line)
                motif = match.group(1)
                if motif not in data:
                    data[motif] = []
            else:
                print('storing data for toyfold3:',motif)
                x = line.strip()
                if float(x) != 0.0: data[motif].append(-log(float(x)))
                #if float(x) != 0.0: data[motif].append(float(x))
    return data


               
main()
