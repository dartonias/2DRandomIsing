import math

tc = math.log(1+2**0.5)/2.

b = [0.1*tc*i for i in range(1,26)]

b = list(set(b))
b = sorted(b,key=lambda x:x)

fout = open('temps.dat','w')
for i in b:
    fout.write('%0.10f\n' % i)
fout.close()
