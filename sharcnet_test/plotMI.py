import numpy as np
import glob

from matplotlib import use
use('agg')

import matplotlib.pyplot as plt


def main():
    files = glob.glob("MI_*.dat")
    t = np.loadtxt('temps.dat')
    for f in files:
        d = np.loadtxt(f)
        L = f.split('_')[1]
        L = L.split('.dat')[0]
        L = int(L.split('L')[1])
        f2 = f.split('MI')[1]
        f2 = 'MIE' + f2
        e = np.loadtxt(f2)
        plt.errorbar(t, d[:,L/2]/L, yerr=e[:,L/2]/L)

if __name__ == "__main__":
    main()
    plt.savefig('test.pdf')
