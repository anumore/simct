#!/usr/bin/env python
import numpy as np
from math import *;
from subprocess import call;
import sys
from glob import glob;
import multiprocessing
import time

## Adds an alpha layer to the final color PNG images (i.e. not visible to the
## naked eye) indicating location of the simulated lens images.
## Aim: To be able to automatically register on which image pixel a user clicked in
## the user interface

call("rm outfits1/CFHTLS_???_????_gri.png",shell=1);

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
        command='./createmask.py %s_o_gri.png blanksims/%s_m_gri.png %s_gri.png'%(subst,subst1,subst);
        call(command,shell=1);

jobs=[];
## Adjust Nproc according to no. of processors on your machine
Nproc=22;
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc))
    jobs.append(p);
    p.start();
