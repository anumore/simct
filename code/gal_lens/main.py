#!/usr/bin/env python
from math import *
import numpy as np
import multiprocessing
import time
import srcprop as srcprop
import srclensprop as slp
from input_qg import *

###########################################################
## INPUTS: 
## Foreground galaxy catalog which has Ra,Dec,Redshift,Magnitudes,
## major-axis,minor-axis,position angle of the ellipticity
##
## PURPOSE:
## Use galaxy catalog to determine which of the galaxies will be potential lenses and
## generate background source magnitudes and redshifts for all sources behind a
## given lensing galaxy.
## Note that background source number density is artificially increased to
## increase the chances of lensing due to the respective foreground galaxy
###########################################################

def worker(num,nproc):
    np.random.seed(num*10+23424);
    srcprop.setlogfile("LOG.%02d.txt"%(num),"GALLENSES.%02d.txt"%(num),"QSOLENSES.%02d.txt"%(num));

    ## Read galaxy lens catalog
    gra,gdec,zd,gu,gg,gr,gi1,gy,gz,majax,minax,ell_pa=np.loadtxt(lenscatalog,usecols=(1,2,5,8,9,10,11,12,13,24,25,22),unpack=True);
    indx,gid,gfld=np.loadtxt(lenscatalog,dtype={'names':('indx','gid','gfld'),'formats':('S3','i8','S13')},usecols=(0,4,3),unpack=True);

    ## Combine i-y band into one
    gi=gi1*1.0;
    for jj in range (gi1.size):
	    if(gi1[jj]<0.):
		    gi[jj]=gy[jj];
    ## Lens ellipticity
    ell=1. - minax/majax;
    qgal=1-ell;

    nqso=0;
    ngal=0;
    iimin=np.int(zd.size*num/nproc);
    iimax=np.int(zd.size*(num+1)/nproc);
    if(num==nproc-1):
        iimax=zd.size;
    for ii in range(iimin,iimax):
        if(ii%10000==0):
            print "Iter",ii,num;
        srcprop.fp1.write("LENSES: %d %d %d %d: \n"%(num,ii,nqso,ngal));
        if(gg[ii]>0 and gr[ii]>0 and zd[ii] >zlens_min and zd[ii] <zlens_max):
            ## Extract Lens Id, Lens Velocity Dispersion, Einstein Radius, Source Magnitude, Source Redshift 
            
            ## For Quasars
            nnew,listmag,listz,rands,vdisp=srcprop.Nsrc_qso(gg[ii],gr[ii],zd[ii],qgal[ii]);
            if(nnew>0 and vdisp>0.):
                for kk in range(nnew):
                    reinst=slp.getreinst(zd[ii],listz[kk],vdisp);
                    ## Keep lenses with Reinst between 1.2 and 5 arcsec only
                if(reinst*206264.8>=reinst_ql and reinst*206264.8 <=reinst_qh):
                    srcprop.fp2.write("%d %f %f %f %f \n"%(ii,vdisp,reinst*206264.8,listmag[kk],listz[kk]));
                    nqso=nqso+1;

            ## For Galaxies
            nnew,listmag,listz,rands,vdisp=srcprop.Nsrc_gal(gg[ii],gr[ii],zd[ii],qgal[ii]);
            if(nnew>0 and vdisp>0.):
                 for kk in range(nnew):
                    reinst=slp.getreinst(zd[ii],listz[kk],vdisp)
                    ## Keep lenses with Reinst between 1.2 and 5 arcsec only
                    if(reinst*206264.8>=reinst_gl and reinst*206264.8 <=reinst_gh):
                        srcprop.fp3.write("%d %f %f %f %f \n"%(ii,vdisp,reinst*206264.8,listmag[kk],listz[kk]));
                        ngal=ngal+1;

    print "Total number of quasar lenses identified by slave:%d in the survey are %d"%(num,nqso)
    print "Total number of galaxy lenses identified by slave:%d in the survey are %d"%(num,ngal)
    srcprop.fp2.close();
    srcprop.fp3.close();
    srcprop.fp1.write("Total number of quasar lenses identified by slave:%d in the survey are %d\n"%(num,nqso));
    tm=time.localtime()
    srcprop.fp1.write("Hour:%d Min:%d Sec:%d"%(tm.tm_hour,tm.tm_min,tm.tm_sec));
    srcprop.fp1.close();

## Run this code faster by specifying Nproc (no. of processors)
jobs=[];
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc))
    jobs.append(p)
    p.start()
print  "#####################################################################";
