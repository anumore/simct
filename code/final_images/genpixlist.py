#!/usr/bin/python
from math import *
import numpy as np
import sys
from StringIO import StringIO

## Extract the x,y pixel position for each lens center
## Copy exposure time from the parent tile to each image cutout 

gid,gra,gdec,fldid,indxno=np.loadtxt('finalpar_srt.txt',dtype={'names': ('gid','gra', 'gdec','fldid','indxno'), 'formats': ('d', 'f', 'f', 'S13', 'S3',)}, usecols=(0,1,2,3,4),unpack=True);

## Create pxlst files with ra,dec list for each tile separately
## Create rgetpix which will run the sky2xy command for all tiles
## Copies exptime from the parent tile to those cutouts which will be finally added to the same parent tile
arr1=np.lexsort((gid,indxno));
gid=gid[arr1];
gra=gra[arr1];
gdec=gdec[arr1];
fldid=fldid[arr1];
indxno=indxno[arr1];

def savet(a,fmt):
    fmt=fmt.split();
    for i in range(len(a)):
        a[i]=np.char.mod(fmt[i],a[i]);

        np.savetxt("inpgal0",np.transpose(a),fmt="%s");

savet([gid,gra,gdec,fldid,indxno],'%d %f %f %s %s');
Ntot=5;
band=['u','g','r','i','z'];

flag=1;
fp1=open('rgetpix','w');
fp2=open('rchd1','w');
for ii in range(gra.size):
    if(flag==1):
        fp=open('imdir/pxlst%s'%(indxno[ii]),'w');
        fp1.write("sky2xy fitsfiles/CFHTLS_W_g_%s_T0007_MEDIAN.fits @imdir/pxlst%s | awk '{print $5,$6}' >> fpxlst\n"%(fldid[ii],indxno[ii]));
	for kk in range(Ntot):
	    fp2.write("cphead fitsfiles/CFHTLS_W_%s_%s_T0007_MEDIAN.fits gout/imoutp_%s_*_%s.fits exptime\n"%(band[kk],fldid[ii],indxno[ii],band[kk]));
        flag=0;
    fp.write("%f %f \n"%(gra[ii],gdec[ii]));
    
    if(ii==(gra.size-1) or indxno[ii]!=indxno[ii+1]):
        fp.close();
        flag=1;

fp1.close();	
fp2.close();	

