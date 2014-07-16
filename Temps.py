bmin = 0.45
bmax = 0.60
N = 30
b = [(bmax-bmin)*i/N + bmin for i in range(N+1)]

#b += [0.1,0.2,0.3,0.4]
b = [0.0,0.1,0.2,0.3,0.4,0.5]

b = list(set(b))
b = sorted(b,key=lambda x:x)

fout = open('temps.dat','w')
for i in b:
    fout.write('%0.6f\n' % i)
fout.close()
