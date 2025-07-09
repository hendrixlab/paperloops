from mpl_toolkits import mplot3d
import numpy as np
import sys, os, re
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from math import log

########
# MAIN #
########

def main():
    usage = "Usage: " + sys.argv[0] + " <loops coords file>"

    if len(sys.argv) != 2:
        print(usage)
        exit(0)

    loop_coords_file = sys.argv[1]

    data = read_loop_coords(loop_coords_file)

    fwd_trace = data[0]
    fwd_helices = data[1]
    rev_trace = data[2]
    rev_helices = data[3]

    fig = plt.figure
    ax = plt.axes(projection='3d')

    for x in fwd_trace:
        start,end,steptype = x
        color = 'blue' if steptype == 'BB' else 'red'
        ax.plot3D([start[0],end[0]],[start[1],end[1]],[start[2],end[2]],color)

    for x in rev_trace:
        start,end,steptype = x
        color = 'blue' if steptype == 'BB' else 'red'
        ax.plot3D([start[0],end[0]],[start[1],end[1]],[start[2],end[2]],color)

    for x in fwd_helices:
        start,end,steptype = x        
        color = 'blue' if steptype == 'BB' else 'red'
        ax.plot3D([start[0],end[0]],[start[1],end[1]],[start[2],end[2]],color)

    for x in rev_helices:
        start,end,steptype = x        
        color = 'blue' if steptype == 'BB' else 'red'
        ax.plot3D([start[0],end[0]],[start[1],end[1]],[start[2],end[2]],color)

    plt.savefig('loop_3D.pdf')
    #plt.show()
        
def read_loop_coords(loop_coords_file):
    data = {}
    with open(loop_coords_file) as F:
        for line in F:
            data = eval(line)            
    return data
                

               
main()
