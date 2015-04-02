"""
Script for submitting jobsto the ASC cluster
for the 'lowa' code, simulating over a range of 
sizes and temperature, and a fixed size of region A

Stephen Inglis, 02.28.2014

"""

import os
import subprocess as sp
import random

def jobNum():
    n = 0
    try:
        fin = open('jobNum','r')
        for i in fin:
            n += 1
        fin.close()
    except IOError:
        pass
    fout = open('jobNum','a')
    fout.write('%d\n' % n)
    fout.close()
    return n

def replace(d):
    r = 'sed'
    for i in d:
        r += " -e's:%s:%s:g'" % (i,d[i])
    return r

def main():
    #pathformat = 'L%03d/b%0.6f'
    pathformat = 'r%03d/R%03d'
    exename = 'Ising.out'
    partemp = 'param_template.dat'
    parname  = 'param.dat'
    subtemp = 'subfile_template'
    subname = 'subfile'
    tempname = 'temps.dat'
    jname = 'Jmat.dat'
    param_dict = {  '__seed':'12345',
                    '__L' : 16,
                    '__ratio' : 0,
                    '__beta':'1.0'}
    sub_dict = {    '__allpath':'cwd_and_exe',
                    '__tarname':'WormTest',
                    '__tardir':'WormTest',
                    '__jobname':'WormTest'}
    jmult = jobNum()
    print 'Batch #%d' % jmult
    jnum = 0
    jmin = 200*jmult
    jmax = 200*(jmult+1)
    #jmin = 0
    #jmax = 10
    disorders = 16
    dis_offset = 0
    #L = 16
    sizein = open('SIZE','r')
    L = int(sizein.readline().split()[1])
    sizein.close()
    for D in range(dis_offset, dis_offset+disorders):
        for R in range(L):
            curr_path = pathformat % (D,R)
            full_path = os.getcwd() + '/' + curr_path
            param_dict['__seed'] = '%d' % random.randint(1,999999)
            param_dict['__beta'] = '%0.8f' % 0.1
            param_dict['__L'] = '%d' % L
            param_dict['__ratio'] = '%d' % R
            sub_dict['__allpath'] = full_path + '/' + exename
            sub_dict['__tarname'] = 'r%03dR%03d.tar' % (D,R)
            sub_dict['__tardir'] = curr_path
            sub_dict['__jobname'] = 'Is2D-%05d' % jnum
            if jnum >= jmin and jnum < jmax:
                try:
                    os.makedirs(curr_path)
                except:
                    pass
                sp.call('cp ../../%s .' % exename,cwd=curr_path,shell=True)
                sp.call(replace(param_dict) + '< ../../%s > %s' % (partemp,parname),cwd=curr_path,shell=True)
                sp.call(replace(sub_dict) + '< ../../%s > %s' % (subtemp,subname),cwd=curr_path,shell=True)
                sp.call('chmod +x %s' % subname,cwd=curr_path,shell=True)
                sp.call('cp ../../%s .' % tempname,cwd=curr_path,shell=True)
                sp.call('cp ../%s .' % jname,cwd=curr_path,shell=True)
                sp.call('touch LOADJ',cwd=curr_path,shell=True)
                sp.call('sqsub --mpp 500M -r 48h -e run.err -o run.out ./subfile',cwd=curr_path,shell=True)
            jnum += 1

if __name__ == "__main__":
    main()
