#!/usr/bin/env python
from math import *
import numpy as np
import sys
from StringIO import StringIO
from subprocess import call
import multiprocessing
import time
import glob
sys.path.append("../");
from input_g import *

###########################################################
## PURPOSE:
## Run lensmodel on gout/*.in with sources as bkg galaxies 
###########################################################

def worker(num,nproc):

    files=glob.glob("ggl*.in")
    iimin=np.int(len(files)*num/nproc);
    iimax=np.int(len(files)*(num+1)/nproc);
    if(num==nproc-1):
	iimax=len(files);
    for ii in range(iimin,iimax):
	print "%s %s > LOG.%s.dat"%(lenscode,files[ii],files[ii]);
	call("%s %s > LOG.%s.dat 2>&1"%(lenscode,files[ii],files[ii]),shell=1);

## Run this code faster by specifying Nproc (no. of processors)
jobs=[];
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc))
    jobs.append(p);
    p.start();
print  "#####################################################################";
