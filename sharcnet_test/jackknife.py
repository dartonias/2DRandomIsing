"""
A script for doing jackknife
"""

import numpy as np

def jackknife(samples, fun):
    # Samples is a list numpy arrays, all of the same length, and assume fun can act on samples of an arbitrary length
    num = len(samples[0])
    jsamp = np.zeros(num)
    for i in xrange(num):
        dat = [np.concatenate((D[:i],D[i+1:])) for D in samples]
        jsamp[i] = fun(dat)
    jjack = jsamp.mean()
    var = 0
    for i in xrange(num):
        var += (jsamp[i] - jjack)**2
    var = var*(num-1.)/num
    return jjack,var**0.5


if __name__ == "__main__":
    print 'Testing jacknife'
    a = np.random.random(100)
    fun = lambda X: X[0].mean()
    print jackknife ([a],fun)
