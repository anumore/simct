#!/usr/bin/env python
import numpy as np
from math import *;
from subprocess import call;
import sys
from glob import glob;
import multiprocessing
import time

## Apply image compression algorithm to final color images

if(len(sys.argv)!=2):
    print "./mkpng_cmb.py dirpath"
    exit(0);
else:    
    dirname=sys.argv[1];

call("rm %s/*png"%(dirname),shell=1);

def worker(num,Nproc):
    fname=glob("outfits1/CFHTLS_*_i.fits");
    iimin=np.int(len(fname)*num/Nproc);
    iimax=np.int(len(fname)*(num+1)/Nproc);
    if(num==Nproc-1):
	iimax=len(fname);

    for ii in range(iimin,iimax):
	fil=fname[ii];
	subst=fil[0:len(fil)-7];
	subst1=fil[9:len(fil)-7];
        command2='pngcrush -d %s %s_gri.png'%(dirname,subst);
        call(command2,shell=1);

jobs=[];
## Adjust Nproc according to no. of processors on your machine
Nproc=22;
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc))
    jobs.append(p);
    p.start();

