#!/usr/bin/env python
from math import *
import numpy as np
import sys
from subprocess import call
import multiprocessing
import time
from input_qg import *

###########################################################
## INPUTS:
## Foreground galaxy catalog (same input as used for main.py), 
## output catalogs from main.py for quasars and galaxies 
##
## PURPOSE:
## Use the lens ids from the all_*.txt and CFHTLS foreground galaxy catalog
## to extract the full lens galaxy catalog with all essential parameters 
###########################################################

print "Reading foreground and background source catalogs and field_ids";
gra,gdec,zd,gu,gg,gr,gi1,gy,gz,majax,minax,ell_pa=np.loadtxt(lenscatalog,usecols=(1,2,5,8,9,10,11,12,13,24,25,22),unpack=True);
indx,gid,gfld=np.loadtxt(lenscatalog,dtype={'names':('indx','gid','gfld'),'formats':('S3','i8','S13')},usecols=(0,4,3),unpack=True);

## Combine i-y band into one
gi=gi1*1.0;
for jj in range (gi1.size):
        if(gi1[jj]<0.):
	    gi[jj]=gy[jj];
## Lens ellipticity
ell=1. - minax/majax;

idxt1,zs1=np.loadtxt("all_qso.txt",usecols=(0,4),unpack=True);
idxt0,zs0=np.loadtxt("all_gal_mod.txt",usecols=(0,4),unpack=True);

print "Initializing arrays";

## Select rows with unique galaxy ids 

## For Quasars
idxtt1,jj= np.unique(idxt1, return_index=True)

idx1=np.arange(idxtt1.size)*0;
sh_str1=np.arange(idxtt1.size)*0.;
sh_pa1=np.arange(idxtt1.size)*0.;

np.random.seed(2892);
for jj in range(idxtt1.size):
    idx1[jj]=int(idxtt1[jj]);
    sh_str1[jj]=np.random.uniform(shear_strength_low,shear_strength_high);
    sh_pa1[jj]=np.random.uniform(shear_pa_low,shear_pa_high);


## For Galaxies
idxtt0,kk= np.unique(idxt0, return_index=True)
idx0=np.arange(idxtt0.size)*0;
sh_str0=np.arange(idxtt0.size)*0.;
sh_pa0=np.arange(idxtt0.size)*0.;

for kk in range(idxtt0.size):
    idx0[kk]=int(idxtt0[kk]);
    sh_str0[kk]=np.random.uniform(shear_strength_low,shear_strength_high);
    sh_pa0[kk]=np.random.uniform(shear_pa_low,shear_pa_high);

print "Found ",jj+kk,"lensing galaxies";

## Extract those galaxies which have one or more background sources within
## the cross-section

## For Quasars
gid1=gid[idx1];
gra1=gra[idx1];
gdec1=gdec[idx1];
gfld1=gfld[idx1];
zd1=zd[idx1];
gu1=gu[idx1];
gg1=gg[idx1];
gr1=gr[idx1];
gi1=gi[idx1];
gz1=gz[idx1];
ell1=ell[idx1];
ell_pa1=ell_pa[idx1];
indxno1=indx[idx1];

## For Galaxies
gid0=gid[idx0];
gra0=gra[idx0];
gdec0=gdec[idx0];
gfld0=gfld[idx0];
zd0=zd[idx0];
gu0=gu[idx0];
gg0=gg[idx0];
gr0=gr[idx0];
gi0=gi[idx0];
gz0=gz[idx0];
ell0=ell[idx0];
ell_pa0=ell_pa[idx0];
indxno0=indx[idx0];


## Save the lens catalogs with background sources as --
## Quasars
np.savetxt('intmpar1.txt',np.transpose((gid1,gra1,gdec1,gfld1,indxno1,zd1,gu1,gg1,gr1,gi1,gz1,ell1,ell_pa1,sh_str1,sh_pa1,idx1)),fmt='%s');
## Galaxies
np.savetxt('intmpar0.txt',np.transpose((gid0,gra0,gdec0,gfld0,indxno0,zd0,gu0,gg0,gr0,gi0,gz0,ell0,ell_pa0,sh_str0,sh_pa0,idx0)),fmt='%s');
call('cat intmpar0.txt intmpar1.txt > intmpar.txt',shell=True);
 
print "Rearranged arrays and kept record of lensing only unique ids of galaxies";
print  "#####################################################################";
sys.stdout.flush();

