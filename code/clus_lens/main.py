#!/usr/bin/env python
from math import *
import numpy as np
import multiprocessing
import time
import srcprop as srcprop
import lum as lum;
from sigma_vd import *
from input_g import *

###########################################################
## INPUTS: 
## Foreground BCG + member catalog which has Ra,Dec,Redshift,Magnitudes,
## major-axis,minor-axis,position angle of the ellipticity
##
## PURPOSE:
## Use galaxy catalog to determine which of the galaxies will be potential lenses and
## generate background source magnitudes and redshifts for all sources behind a
## given lensing galaxy.
## Note that background source number density is artificially increased to
## increase the chances of lensing due to the respective foreground galaxy
###########################################################

## Read galaxy lens catalog
gra,gdec,zd,gu,gg,gr,gi,gz,majax,minax,ell_pa,mid=np.loadtxt(lenscatalog,usecols=(1,2,5,8,9,10,11,13,24,25,22,31),unpack=True);
indx,gid,gfld=np.loadtxt(lenscatalog,dtype={'names':('indx','gid','gfld'),'formats':('S3','i8','S13')},usecols=(0,4,3),unpack=True);

## Lens velocity dispersion
sig=np.arange(zd.size)*0.;
for jj in range (gra.size):
    sig[jj]=getsigma(gg[jj],gr[jj],zd[jj]);

## Lens ellipticity
ell=1. - minax/majax;

## Requires catalog file to be sorted in ascending order
nwmid,bcgid= np.unique(mid, return_index=True)

def worker(num,nproc):
    np.random.seed(num*10+23424);
    srcprop.setlogfile("LOG.%02d.txt"%(num),"GALLENSES.%02d.txt"%(num),"LENCAT.%02d.txt"%(num));
    
    ngal=0;
    iimin=np.int(bcgid.size*num/nproc);
    iimax=np.int(bcgid.size*(num+1)/nproc);
    if(num==nproc-1):
	iimax=bcgid.size;
    

    for ii in range(iimin,iimax):
	#print "iter",ii,num;
        
	srcprop.fp1.write("LENSES: %d %d %d: "%(num,ii,ngal));
         
	## Calculate kappa excluding the source redshift dependent factor
	## hinv Mpc, hinv Msun/(hinv Mpc)^3
        rscale,rhos=lum.getRhos_Rs_i(gi[bcgid[ii]],zd[bcgid[ii]]);
        ## Msun/(hinv Mpc^2)
	kfac=rscale*rhos;
        Dd=cc.Daofz(zd[bcgid[ii]]);
        ## Msun/(hinv Mpc^2)
	sigcrit=cspd**2/(4*pi*gee*Dd);  
	##dimensionless
	kappa1=kfac/sigcrit;
	
	minv=bcgid[ii];
	if(ii==nwmid.size-1):
	    maxv=mid.size;
	else:
	    maxv=bcgid[ii+1];
        
	## BCG+member positions in pixels
	gposx=(gra[bcgid[ii]]-gra[minv:maxv])*3600/pixsc;
	gposy=(-gdec[bcgid[ii]]+gdec[minv:maxv])*3600/pixsc;
	rscale=(rscale/Dd)*(206264.8/pixsc);
        
	## Save a catalog with BCG+member parameters
	for tt in range(minv,maxv):
	    ll=tt-minv;
	    srcprop.fp3.write("%d %d %f %f %f %f %f %f %f %f\n"%(mid[tt],gid[tt],kappa1,rscale,gposx[ll],gposy[ll],zd[tt],sig[tt],ell[tt],ell_pa[tt]));
        
	## Calculate and save the magnitudes and redshifts of multiple sources per lens 
        nnew,listmag,listz,rands=srcprop.Nsrc_gal(kappa1,rscale,sig[minv:maxv],gposx,gposy,ell[minv:maxv],ell_pa[minv:maxv],zd[minv],ii,num);
        if(nnew>0):
            for kk in range(nnew):
                srcprop.fp2.write("%d %d %d %f %d %f %f \n"%(int(nwmid[ii]),gid[bcgid[ii]],bcgid[ii],zd[bcgid[ii]],nnew,listmag[kk],listz[kk]));
                ngal=ngal+1;

    srcprop.fp1.write("Total number of galaxy lenses identified by slave:%d in the survey are %d\n"%(num,ngal));
    print "Total number of background lensed galaxies identified by slave:%d in the survey are %d"%(num,ngal)
    tm=time.localtime()
    srcprop.fp1.write("Hour:%d Min:%d Sec:%d"%(tm.tm_hour,tm.tm_min,tm.tm_sec));
    srcprop.fp1.close();
    srcprop.fp2.close();
    srcprop.fp3.close();

## Run this code faster by specifying Nproc (no. of processors)
jobs=[];
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc))
    jobs.append(p)
    p.start()
print  "#####################################################################";
