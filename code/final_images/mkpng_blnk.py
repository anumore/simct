#!/usr/bin/env python
import numpy as np
from math import *;
from subprocess import call;
import sys
from glob import glob;
import multiprocessing
import time

## Creates color PNG images of the simulated lensed images only i.e. no real
## survey images

## Needs compose.py from HumVI package

call("rm blanksims/*png",shell=1);

def worker(num,Nproc):
    fname=glob("blanksims/CFHTLS_*_i.fits");

    iimin=np.int(len(fname)*num/Nproc);
    iimax=np.int(len(fname)*(num+1)/Nproc);
    if(num==Nproc-1):
	iimax=len(fname);

    for ii in range(iimin,iimax):
	fil=fname[ii];
	subst=fil[0:len(fil)-7];
	command='./compose.py -s 0.4,0.6,1.7 -z 0.0 -p 1.0,0.09 -m -1.0  -o %s_m_gri.png %s_i.fits %s_r.fits %s_g.fits'%(subst,subst,subst,subst);
	#print command
        call(command,shell=1);


jobs=[];
## Adjust Nproc according to no. of processors on your machine
Nproc=22;
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc))
    jobs.append(p);
    p.start();

