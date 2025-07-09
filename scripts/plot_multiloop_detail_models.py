import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys, os
from scipy.optimize import curve_fit
from matplotlib.colors import LinearSegmentedColormap
import random
from sklearn.metrics import r2_score

def main():
    usage = "Usage: " + sys.argv[0] + " <asymmetry values>"
    if len(sys.argv) != 2:
        print(usage)
        exit()

    data_file = sys.argv[1]    
    data = read_data(data_file)
    
    for n in data:
        for L in data[n]:
            if int(n) > 2 and int(L) > 0:
                #print("examining:",n,L)
                labels,asyms,balances,skews,lengths,dGs = zip(*data[n][L])
                data1,data2 = get_coeff(labels)
                asyms = np.array(asyms)
                balances = np.array(balances)
                skews = np.array(skews)
                lengths = np.array(lengths)
                x_values = (asyms,balances,data1,data2)
                y_values = np.array(dGs)
                if len(y_values) >= 5:
                    #print(y_values,x_values)
                    popt,pcov = curve_fit(f,x_values,y_values)
                    y_pred = f(x_values,*popt)
                    R2 = r2_score(y_pred,y_values)
                    (A,B,C,D,E) = popt
                    print(n,L,R2,A,B,C,D,E)

#############
# FUNCTIONS #
#############

def get_coeff(labels):
    data1 = []
    data2 = []
    for label in labels:
        good = True
        lengths = [float(this_L) for this_L in label.split('-') ]
        fc = np.fft.fft(lengths)
        if sum(fc) == 0:
            good = False
        if good:
            #print(fc)
            fc2 = np.abs(fc)**2
            #print(fc2)
            fc2 = fc2 / np.sum(fc2)
            P1 = fc2[int((len(fc2)-1)/2)]
            P2 = fc2[len(fc2)-1]
            data1.append(float(P1))
            data2.append(float(P2))
        else:
            data1.append(0.0)
            data2.append(0.0)

    return np.array(data1),np.array(data2)
    
def f(x,A,B,C,D,E):
    asym,H,P1,P2 = x
    return  A*asym + B*H + C*P1 + D*P2 + E

def read_data(data_file):
    F = open(data_file)
    data = {}
    for line in F:
        n,L,label,asym,balance,skew,dG = line.strip().split('\t')
        if n not in data:
            data[n] = {}
        if L not in data[n]:
            data[n][L] = []
        if not np.isnan(float(dG)):
            data[n][L].append((label,float(asym),float(balance),float(skew),int(L),float(dG)))
    return data

main()
