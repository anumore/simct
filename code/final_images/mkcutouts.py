#!/opt/local/bin/python2.7
from math import *
import numpy as np
import sys
from StringIO import StringIO
from subprocess import call

## Generate a file which upon running will make cutouts of given size and offsets for each tile in the input file
## Also, creates cutoutlist

## width of each cutout
wd=440;

## pixel offset for between neighbouring cutouts
## (since we would like to  have some overlapping region; hence dx!=wd)
dx=386;

## Total number of cutouts is Ntot**2
Ntotpix=19354;
Ntot=Ntotpix/dx;

if(len(sys.argv)!=5):
    print "./mkcutouts.py dirname width offset ntotpix \n dirname: path of the output
    directory, width of each cutout in pixels, spatial offset in pixels between
    neighbouring cutouts, total no. of pixels in a given parent tile";
    sys.exit(0);
else:
    dirname=sys.argv[0];
    wd=int(sys.argv[1]);
    dx=int(sys.argv[2]);
    Ntotpix=int(sys.argv[3]);

## Read the tile IDs and names
indxno,fldid=np.loadtxt('fieldid_sort',dtype={'names':('indxno','fldid'),'formats':('S9','S13')},usecols=(0,1),unpack=True);

band=['u','g','r','i','z'];
Ntotb=len(band);

## A combined executable file for all tiles
fp1=open('imdir/rcoutall','w');
fp1.write('#!/bin/bash \n');

## A catalog which has: "xmin,xmax,ymin,ymax,cutout_no." wrt the parent tile 
fp2=open('cutoutlist','w');

for kk in range(fldid.size):
    ## Executable file for each tile
    fp=open('imdir/rcutout_'+indxno[kk],'w');
    fp.write('#!/bin/bash \n');
    for ii in range(Ntot):
	for jj in range(Ntot):
	    for ll in range(Ntotb):
		fp.write('fitscopy.out fitsfiles/CFHTLS_W_%s_%s_T0007_MEDIAN.fits[%d:%d,%d:%d] %s/CFHTLS_%s_%04d_%s.fits \n' %(band[ll],fldid[kk],(ii*dx)+1,(ii*dx)+wd, (jj*dx)+1,(jj*dx)+wd,dirname,indxno[kk],jj+1+(Ntot*ii),band[ll]));
		fp2.write('%d %d %d %d %04d \n'%((ii*dx)+1,(ii*dx)+wd, (jj*dx)+1,(jj*dx)+wd,jj+1+(Ntot*ii)));
    fp.close(); 
    fp1.write('imdir/rcutout_%s \n'%(indxno[kk]));
fp1.close(); 
fp2.close(); 

call("chmod +x imdir/rcutout_* imdir/rcoutall",shell=1);
