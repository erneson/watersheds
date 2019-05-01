#!/usr/bin/env python3

import sys
import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import random as rand

plt.rc('patch',linewidth=2)
plt.rc('axes', linewidth=2, labelpad=10)
plt.rc('xtick.minor', size=4, width=2)
plt.rc('xtick.major', size=8, width=2, pad=8)
plt.rc('ytick.minor', size=4, width=2)
plt.rc('ytick.major', size=8, width=2, pad=8)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern', size=30)

name=sys.argv[1]
hcutoff=sys.argv[2]

f=open('sigma_'+name+'_h'+hcutoff+'.dat','r')
line=f.readline()
field=line.split()
lx=int(field[0])
ly=int(field[1])
n=lx*ly

ncol=lx
nrow=ly
sig=-2*np.ones((nrow,ncol),dtype=int)
for line in f:
    field=line.split()
    
    i=int(field[0])
    j=int(field[1])
    s=int(field[2])
    
    sig[j,i]=s
f.close()

mydict={}
for j in range(ly):
    for i in range(lx):
        s=sig[j,i]
        
        if s in mydict:
            mydict[s].append((j,i))
        else:
            mydict[s]=[(j,i)]

ratio=float(ly)/float(lx)

xbox=10
ybox=int(xbox*ratio)
fig,ax=plt.subplots(figsize=(xbox,ybox),dpi=300)

N=len(mydict)
keys=sorted(list(mydict.keys()))

t=0
for k in keys:
    mylist=mydict[k]
    
    for mytuple in mylist:
        j=mytuple[0]
        i=mytuple[1]
        
        sig[j,i]=t
    
    t=t+1

rand.seed(1231)
if(keys[1]==-1):
    mycolors = [(0.9,0.9,0.9)]+[(0.0,0.0,0.0)]+[(0.2+0.5*rand.random(),0.2+0.5*rand.random(),0.2+0.5*rand.random()) for i in range(N-2)]
else:
    mycolors = [(0.9,0.9,0.9)]+[(0.2+0.5*rand.random(),0.2+0.5*rand.random(),0.2+0.5*rand.random()) for i in range(N-1)]

mymap = mpl.colors.LinearSegmentedColormap.from_list('mymap', mycolors, N=N)

myplot=plt.imshow(sig,cmap=mymap,interpolation='none')
myplot.axes.get_xaxis().set_visible(False)
myplot.axes.get_yaxis().set_visible(False)
plt.axis('off')

plt.savefig('sigma_'+name+'_h'+hcutoff+'.pdf', bbox_inches='tight')
