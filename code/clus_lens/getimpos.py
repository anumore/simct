#!/usr/bin/env python
from math import *
import numpy as np
from subprocess import call
import linecache
import string
from input_g import *

###########################################################
## INPUTS:
## Catalog of lensed sources which has filename, source positions and lens id
## and another catalog with source positions, no. of lensed images, row number
## and lens id
##
## PURPOSE:
## Extract those lenses which have Reinst > Reinst_min where Reinst_min=2.0
## arcsecond
###########################################################

## Read source info and lens id
fname,gid=np.loadtxt('srcinf_inp',dtype={'names':('fname','gid'),'formats':('S40','d11')},usecols=(0,5),unpack=1);
sx,sy=np.loadtxt('srcinf_inp',usecols=(1,2),unpack=1);

sx1,sy1,imgno,rwn,gid1=np.loadtxt('srcinfoall',usecols=(0,1,2,3,4),unpack=1);

ngid1,uid=np.unique(gid1,return_index=True);

fidx=np.arange(gid.size)*0;

## Extract the lensed image positions and magnifications corresponding to each
## source and calculate the Reinst
call('rm imglist',shell=1);
nn=0;
for ii in range(gid.size):
    for jj in range(ngid1.size):
	if(gid[ii]-ngid1[jj]==0):
            kk=uid[jj];
	    cond=True;
	    while(cond):
		## Store the image x,y-positions, magnifications, lensing galaxy
		## id and source x,y-positions for the full sample
		if(sx[ii]-sx1[kk]==0 and sy[ii]-sy1[kk]==0):
		    call("awk \'{if(NR>='%d' && NR<'%d') print $1,$2,$3,'%d','%f','%f'"%(rwn[kk],rwn[kk]+imgno[kk],gid[ii],sx[ii],sy[ii])+"}\' %s >> imglist" %(fname[ii]),shell=1);
		    call("awk \'{if(NR>='%d' && NR<'%d') print $1,$2,$3,'%d'"%(rwn[kk],rwn[kk]+imgno[kk],gid[ii])+"}\' %s > imginfo_t" %(fname[ii]),shell=1);
		
		    ## Store the image x,y-positions and magnifications for a
		    ## lens
		    imx,imy,immag=np.loadtxt('imginfo_t',usecols=(0,1,2),unpack=1);
                    
		    idxn=np.argmin(imx);
		    idxx=np.argmax(imx);
		    sep1=np.sqrt( (imx[idxn]-imx[idxx])**2+ (imy[idxn]-imy[idxx])**2 );
		    
		    idyn=np.argmin(imy);
		    idyx=np.argmax(imy);
		    sep2=np.sqrt( (imx[idyn]-imx[idyx])**2+ (imy[idyn]-imy[idyx])**2 );

		    Reinst=pixsc*((sep1+sep2)/2.)/2.;
		    ## Reinst_min in arcsec
		    if(Reinst>Reinst_min):
			fidx[nn]=ii;
			nn=nn+1;
		    cond=False;
		else:
		    kk=kk+1;


## Store ids of those lenses which have Reinst>Reinst_min
fidx=fidx[0:nn];
nnmax=nn;
print "Found",nnmax,"cluster-scale lenses with Reinst >",Reinst_min,"arcsec";

## Extract the full catalog such that lenses have Reinst>Reinst_min
cur=0
nn=0
fp1=open('finalpar0.txt','r');
fp2=open('finpar_g2.txt','w');
for line in fp1:
    line=line[:-1];
    if(nn<nnmax and cur==fidx[nn]):
	fp2.write("%s \n"%(line));
	nn=nn+1;
    cur=cur+1;
fp2.close();

## Extract randomly any Ncl_lens lenses out of this sample
## E.g. 630 is the number of lenses we decided to simulate in CFHTLS-W1 field
np.random.seed(7238);
fp2=open('finalpar.txt','w');
filename='finpar_g2.txt';
rni=np.random.permutation(fidx.size)+1;
written=0;
for kk in range(rni.size):
    line=linecache.getline(filename,rni[kk])[:-1];
    values=string.split(line,sep=None);
    magi=float(values[9]);
    if(magi<magi_upperlim):
	continue;
    fp2.write("%s \n"%(line));
    written=written+1;
    if(written==Ncl_lens):
	break;

print "Finally, extracting a random set of",Ncl_lens,"cluster-scale lenses out of",nnmax;
print  "#####################################################################";
