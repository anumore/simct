#!/usr/bin/env python
from math import *
import numpy as np
import sys
from StringIO import StringIO
from subprocess import call
import multiprocessing
import time
import glob
import os;

print "Deleting blanksims/*fits outfits/*fits outfits1/*fits imdir/LOG*.dat";
sys.stdout.flush();
call("rm ../blanksims/CF*fits",shell=1); 
call("rm ../outfits/CF*fits",shell=1); 
call("rm ../outfits1/CF*fits",shell=1); 
call("rm LOG*.dat",shell=1); 

def worker(num,nproc):

    files=glob.glob("rmgcpy_*")
    files1=glob.glob("rmrg_*")
    files2=glob.glob("rmrgb_*")
    files3=glob.glob("rimst_*")
    iimin=np.int(len(files)*num/nproc);
    iimax=np.int(len(files)*(num+1)/nproc);
    if(num==nproc-1):
        iimax=len(files);
    for ii in range(iimin,iimax):
        print "Running copy %s"%(files[ii]);
        sys.stdout.flush();
        call("./%s 2>&1"%(files[ii]),shell=1);
   
    fid,gid=np.loadtxt("idpxlst",dtype={'names':('fid','gid'),'formats':('S8','d')},usecols=([0,1]),unpack=1);
    done=False;
    cntloop=0;
    while(not done):
	time.sleep(3);
	done=True;
	for ii in range(fid.size):
	    for col in ['u','g','r','i','z']:
		ff="../outfits/CFHTLS_%s_%s_t.fits"%(fid[ii],col);
		if(not os.path.isfile(ff)):
		    if(cntloop%1000==0):
			print "File",ff," does not exist"
		    done=False;
	cntloop=cntloop+1;
	if(cntloop>10000):
	    print "Some files missing in outfits1, exiting";
	    exit(0);
        
    for ii in range(iimin,iimax):
        print "Running merge %s"%(files[ii]);
        sys.stdout.flush();
        call("./%s > LOG.%s.dat 2>&1 "%(files1[ii],files1[ii]),shell=1);
        
    for ii in range(iimin,iimax):
        print "Running mergeblank%s"%(files[ii]);
        sys.stdout.flush();
        call("./%s > LOG.%s.dat 2>&1"%(files2[ii],files2[ii]),shell=1);

    for ii in range(iimin,iimax):
        print "Running imstat %s"%(files[ii]);
        sys.stdout.flush();
        call("./%s > LOG.%s.dat 2>&1"%(files3[ii],files3[ii]),shell=1);

jobs=[];
Nproc=22;
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc))
    jobs.append(p);
    p.start();

