import subprocess as sp
import os

def main():
    sizein = open('SIZE','r')
    L = int(sizein.readline().split()[1])
    sizein.close()
    folder_f = 'L%d' % L
    folder_r = 'L%dr' % L
    try:
        os.makedirs(folder_f)
        os.makedirs(folder_r)
    except:
        pass
    sp.call('python genJs.py',shell=True)
    sp.call('tar -cf jmat.tar r*/Jmat.dat',shell=True)
    sp.call('rm -rf r*/',shell=True)
    files = ['jmat.tar', 'Ising.out', 'param_template.dat', 'SIZE', 'subfile_template', 'submit_sharcnet.py', 'tar_stuff', 'temps.dat']
    for f in files:
        sp.call('cp ../%s .' % f,cwd=folder_f,shell=True)
        sp.call('cp ../%s .' % f,cwd=folder_r,shell=True)
    sp.call('tar -xf jmat.tar',cwd=folder_f,shell=True)
    sp.call('tar -xf jmat.tar',cwd=folder_r,shell=True)
    sp.call('rm jmat.tar',cwd=folder_f,shell=True)
    sp.call('rm jmat.tar',cwd=folder_r,shell=True)

if __name__ == "__main__":
    main()
