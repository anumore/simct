#!/usr/bin/env python
from math import *
import numpy as np

###########################################################
## INPUTS:
## output source catalogs from main.py
##
## PURPOSE:
## uniqidx has lens id for source1 catalog
## all_xxx.txt is the source2 catalog from which lens ids common to source1
## catalog are removed  
###########################################################

## Uncomment following 3 lines, if you want to remove the common lensing galaxy records from the
## bkg gal catalog
idxt1,sigm1,reinst1,smii1,zs1=np.loadtxt("all_gal.txt",unpack=True);
idg=np.loadtxt('uniqidq',unpack=True);
fp=open('all_gal_mod.txt','w');

## Uncomment following 3 lines, if you want to remove the common lensing galaxy records from the
## bkg qso catalog

#idxt1,sigm1,reinst1,smii1,zs1=np.loadtxt("all_qso.txt",unpack=True);
#idg=np.loadtxt('uniqidg',unpack=True);
#fp=open('all_qso_mod.txt','w');

ii=0;
jj=0;
while jj < idg.size and ii<idxt1.size:
	if(idxt1[ii]<=idg[jj]):
	    if(idxt1[ii]==idg[jj]):
		ii=ii+1;
	    else:
		fp.write(" %d %f %f %f %f \n"%(idxt1[ii],sigm1[ii],reinst1[ii],smii1[ii],zs1[ii]));
		ii=ii+1;
	else:
	    jj=jj+1;
	
