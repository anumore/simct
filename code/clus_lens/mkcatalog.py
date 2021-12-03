#!/usr/bin/env python
from math import *
import numpy as np
import sys
from subprocess import call
import multiprocessing
import time
from input_g import *

###########################################################
## INPUTS:
## Foreground cluster catalog (same input as used for main.py), 
## output catalogs from main.py for bkg galaxies 
##
## PURPOSE:
## Use the lens ids from the all_*.txt and CFHTLS foreground galaxy catalog
## to extract the full lens galaxy catalog with all essential parameters 
###########################################################


## Use the lens ids from the all_*.txt and CFHTLS foreground galaxy catalog
## to extract the full lens galaxy catalog with all essential parameters 
print "Reading foreground and background source catalogs and field_ids";
gra,gdec,zd,gu,gg,gr,gi,gz,majax,minax,ell_pa,mid=np.loadtxt(lenscatalog,usecols=(1,2,5,8,9,10,11,13,24,25,22,31),unpack=True);

indx,gid,gfld=np.loadtxt(lenscatalog,dtype={'names':('indx','gid','gfld'),'formats':('S3','i8','S13')},usecols=(0,4,3),unpack=True);

ell=1. - minax/majax;

np.random.seed(2892);
print "Initializing arrays";

ngid,ouid=np.loadtxt("all_gal.txt",dtype={'names':('ngid','ouid'),'formats':('i8','i8')},usecols=(1,2),unpack=True);


nwid,nouid=np.unique(ouid,return_index=True);
print "Found ",nouid.size,"lenses";

## Select rows with unique galaxy ids 
gid0=gid[nwid];
gra0=gra[nwid];
gdec0=gdec[nwid];
gfld0=gfld[nwid];
zd0=zd[nwid];
gu0=gu[nwid];
gg0=gg[nwid];
gr0=gr[nwid];
gi0=gi[nwid];
gz0=gz[nwid];
ell0=ell[nwid];
ell_pa0=ell_pa[nwid];
indxno0=indx[nwid];

np.savetxt('intmpar.txt',np.transpose((gid0,gra0,gdec0,gfld0,indxno0,zd0,gu0,gg0,gr0,gi0,gz0,ell0,ell_pa0)),fmt='%s');
 
print "Rearranged arrays and kept record of lensing only unique ids of galaxies";
print "###############################################################";
sys.stdout.flush();

