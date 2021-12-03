#!/usr/bin/env python
from math import *
import numpy as np
import sys
from StringIO import StringIO
from subprocess import call
from math import *
import multiprocessing
import time
import glob 

def worker(num,nproc):
    files=glob.glob("gout/imoutp*.fits")
    iimin=np.int(len(files)*num/nproc);
    iimax=np.int(len(files)*(num+1)/nproc);
    if(num==nproc-1):
        iimax=len(files);
    for ii in range(iimin,iimax):	
	call("./poisson.py %s 2>&1"%(files[ii]),shell=1);

jobs=[];
## Adjust Nproc according to no. of processors on your machine
Nproc=22;
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc))
    jobs.append(p);
    p.start();
