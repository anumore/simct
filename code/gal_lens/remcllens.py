#!/usr/bin/env python
from math import *
from subprocess import check_output
import numpy as np
import sys
import string
import linecache

## Read the cluster lens catalog
clsid=np.loadtxt(cluslenscatalog,usecols=([0]),unpack=True);

## For Quasars
fp=open('finalpar1.txt','r');
fp2=open('finalpar1_clsb.txt','w');

## Extract those lenses which do not match with the lenses in the cluster lens catalog
jj=0;
for line in fp:
    ## Read quasar id
    line=line[:-1];
    values=string.split(line,sep=None);
    gid=int(values[0]);
    
    ## Advance jj until you have its index>=gid
    while(gid>clsid[jj] and jj<clsid.size-1):
	jj=jj+1;
 
    if(gid!=clsid[jj]):
	fp2.write("%s \n"%(line));

fp.close();
fp2.close();

## Extract randomly any N_lens lenses out of this sample
## E.g. 630 is the number of lenses we decided to simulate in CFHTLS-W1 field
fp2=open('finalpar_q.txt','w');
filename='finalpar1_clsb.txt';
fp=open(filename,"r");
qsolim=int(check_output("wc -l finalpar1_clsb.txt",shell=1).split()[0]);
rni=np.random.permutation(qsolim);
idx=np.lexsort([rni]);
cnt=0;
## Note N_lens must be < qsolim
for line in fp:
    line=line[:-1];
    if(idx[cnt]<N_lens):
	fp2.write("%s \n"%(line));
    cnt=cnt+1;

fp.close();
fp2.close();

## For Galaxies
fp=open('finalpar0.txt','r');
fp2=open('finalpar0_clsb.txt','w');

## Extract those lenses which do not match with the lenses in the cluster lens catalog
jj=0;
for line in fp:
    ## Read galaxy id
    line=line[:-1];
    values=string.split(line,sep=None);
    gid=int(values[0]);
    
    ## Advance jj until you have its index>=gid
    while(gid>clsid[jj] and jj<clsid.size-1):
	jj=jj+1;

    if(gid!=clsid[jj]):
	fp2.write("%s \n"%(line));

fp.close();
fp2.close();

## Extract randomly any N_lens lenses out of this sample
## E.g. 630 is the number of lenses we decided to simulate in CFHTLS-W1 field
fp2=open('finalpar_g.txt','w');
filename='finalpar0_clsb.txt';
fp=open(filename,"r");
gallim=int(check_output("wc -l finalpar0_clsb.txt",shell=1).split()[0]);
rni=np.random.permutation(gallim);
idx=np.lexsort([rni]);
cnt=0;
## Note N_lens must be < gallim
for line in fp:
    line=line[:-1];
    if(idx[cnt]<N_lens):
	fp2.write("%s \n"%(line));
    cnt=cnt+1;

fp2.close();
fp.close();

