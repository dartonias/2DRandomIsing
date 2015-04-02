"""
Analysis script for 2D Ising model
"""

import tarfile
import numpy as np
from jackknife import jackknife as jk

def loadTemps():
    r = []
    fin = open('temps.dat')
    for i in fin:
        r.append(float(i))
    return r

def analyze(L=8, ex=''):
    filef = "L%d%s/data.tar"
    filer = "L%d%sr/data.tar"
    fmt = "r%03d/R%03d/obs_ratio_%01.6f"
    betas = loadTemps()
    disorders = range(128)
    sizes = range(0,L)
    
    dataf = tarfile.open(filef % (L,ex),'r')
    t_entf = [[[0. for i in xrange(L+1)] for z in disorders] for zz in betas]
    for Bn,BB in enumerate(betas):
        for DD in disorders:
            for SS in sizes:
                f = np.loadtxt(dataf.extractfile(fmt % (DD,SS,BB)))
                val = f.mean()
                t_entf[Bn][DD][SS+1] += -1*np.log(val) + t_entf[Bn][DD][SS]
    dataf.close()

    datar = tarfile.open(filer % (L,ex),'r')
    t_entr = [[[0. for i in xrange(L+1)] for z in disorders] for zz in betas]
    for Bn,BB in enumerate(betas):
        for DD in disorders:
            for SS in sizes:
                f = np.loadtxt(datar.extractfile(fmt % (DD,SS,BB)))
                val = f.mean()
                t_entr[Bn][DD][SS+1] += -1*np.log(val) + t_entr[Bn][DD][SS]
    datar.close()

    t_mi = [[[0. for i in xrange(L+1)] for z in disorders] for zz in betas]
    for Bn,BB in enumerate(betas):
        for DD in disorders:
            for SS in sizes:
                t_mi[Bn][DD][SS+1] = 0.5*(
                (t_entf[Bn][DD][SS+1] + t_entr[Bn][DD][SS+1]) +
                (t_entf[Bn][DD][-1*SS-2] + t_entr[Bn][DD][-1*SS-2]) -
                (t_entf[Bn][DD][-1] + t_entr[Bn][DD][-1])
                )

    t_mi = np.array(t_mi)
    final = [[0. for i in xrange(L+1)] for zz in betas]
    finalE = [[0. for i in xrange(L+1)] for zz in betas]
    func = lambda X: X[0].mean()
    for Bn,BB in enumerate(betas):
        for SS in sizes:
            val, err = jk([t_mi[Bn,:,SS+1]],func)
            final[Bn][SS+1] = val
            finalE[Bn][SS+1] = err
    np.savetxt('MI_L%d%s.dat' % (L,ex), np.array(final))
    np.savetxt('MIE_L%d%s.dat' % (L,ex), np.array(finalE))

if __name__ == "__main__":
    analyze(8,'g')
    analyze(12,'g')
