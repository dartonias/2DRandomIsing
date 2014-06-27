"""
Script for submitting jobsto the ASC cluster
for the 'lowa' code, simulating over a range of 
sizes and temperature, and a fixed size of region A

Stephen Inglis, 02.28.2014

"""

import os
import subprocess as sp
import random

def replace(d):
    r = 'sed'
    for i in d:
        r += " -e's:%s:%s:g'" % (i,d[i])
    return r

def main():
    pathformat = 'r%03d/R%03d'
    exename = 'Ising.out'
    partemp = 'param_template.dat'
    parname  = 'param.dat'
    subtemp = 'subfile_template'
    subname = 'subfile'
    Jfile = 'Jmat.dat'
    param_dict = {  '__seed':'12345',
                    '__regionA':'0'}
    sub_dict = {    '__allpath':'cwd_and_exe',
                    '__tarname':'WormTest',
                    '__tardir':'WormTest',
                    '__jobname':'WormTest'}
    jmult = 0
    jnum = 0
    jmin = 128*jmult
    jmax = 128*(jmult+1)
    #jmin = 1
    #jmax = 128
    for r in xrange(1):
        for R in xrange(64,256,4):
            curr_path = pathformat % (r,R)
            full_path = os.getcwd() + '/' + curr_path
            param_dict['__seed'] = '%d' % random.randint(1,999999)
            param_dict['__regionA'] = '%d' % R
            sub_dict['__allpath'] = full_path + '/' + exename
            sub_dict['__tarname'] = 'r%03dR%03d.tar' % (r,R)
            sub_dict['__tardir'] = curr_path
            sub_dict['__jobname'] = 'SK-%04d' % jnum
            if jnum >= jmin and jnum < jmax:
                try:
                    os.makedirs(curr_path)
                except:
                    pass
                sp.call('cp ../../%s .' % exename,cwd=curr_path,shell=True)
                sp.call('cp ../%s .' % Jfile,cwd=curr_path,shell=True)
                sp.call(replace(param_dict) + '< ../../%s > %s' % (partemp,parname),cwd=curr_path,shell=True)
                sp.call(replace(sub_dict) + '< ../../%s > %s' % (subtemp,subname),cwd=curr_path,shell=True)
                sp.call('qsub subfile',cwd=curr_path,shell=True)
            jnum += 1

def main2():
    betas = [0.1,0.2]
    pathformat = 'b%0.6f'
    exename = 'Ising.out'
    partemp = 'param_template.dat'
    parname  = 'param.dat'
    subtemp = 'subfile_template'
    subname = 'subfile'
    param_dict = {  '__seed':'12345',
                    '__regionA':'0'}
    sub_dict = {    '__allpath':'cwd_and_exe',
                    '__tarname':'WormTest',
                    '__tardir':'WormTest',
                    '__jobname':'WormTest'}
    jmult = 0
    jnum = 0
    jmin = 128*jmult
    jmax = 128*(jmult+1)
    for b in betas:
        curr_path = pathformat % b
        full_path = os.getcwd() + '/' + curr_path
        param_dict['__seed'] = '%d' % random.randint(1,999999)
        sub_dict['__allpath'] = full_path + '/' + exename
        sub_dict['__tarname'] = 'b%0.6f.tar' % (r,R)
        sub_dict['__tardir'] = curr_path
        sub_dict['__jobname'] = '2DIs-%04d' % jnum
        if jnum >= jmin and jnum < jmax:
            try:
                os.makedirs(curr_path)
            except:
                pass
            sp.call('cp ../%s .' % exename,cwd=curr_path,shell=True)
            sp.call(replace(param_dict) + '< ../%s > %s' % (partemp,parname),cwd=curr_path,shell=True)
            sp.call(replace(sub_dict) + '< ../%s > %s' % (subtemp,subname),cwd=curr_path,shell=True)
            sp.call('qsub subfile',cwd=curr_path,shell=True)
        jnum += 1

if __name__ == "__main__":
    main()
