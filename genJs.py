"""
Script for submitting jobsto the ASC cluster
for the 'lowa' code, simulating over a range of 
sizes and temperature, and a fixed size of region A

Stephen Inglis, 02.28.2014

"""

import os
import subprocess as sp
import numpy as np

def main():
    #pathformat = 'L%03d/b%0.6f'
    pathformat = 'r%03d'
    disorders = 128
    L = 16
    for D in range(0,disorders+0):
        curr_path = pathformat % D
        try:
            os.makedirs(curr_path)
        except OSError:
            pass
        fout = open(curr_path + '/Jmat.dat','w')
        for i in xrange(L**2 * 2):
            fout.write("%1.8f " % np.random.normal())
        fout.close()

if __name__ == "__main__":
    main()
