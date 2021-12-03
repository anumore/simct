#!/usr/bin/env python
from math import *
from subprocess import call
import multiprocessing
import time
import numpy as np
import sys
import gensrcpos as scp
import genimg as gn
from input_qg import *

###########################################################
## INPUTS: 
## Foreground galaxy catalog with lens model properties as output by
## mkcatalog.py and Background quasar catalog which has Redshift and Magnitude 
##
## PURPOSE:
## Use the lens properties to choose a single source per lens, calculate 
## source position based on the limits on fluxes of the lensed images, assign
## colors and generate final lensed quasar images 
##
###########################################################

## Quasars as background sources

def worker(num,nproc):
    ## Read the lens galaxy catalog, real bkg quasar color catalog and 
    ## background quasar catalog generated for each potential lens 
    gid1,gra1,gdec1,zd1,gu1,gg1,gr1,gi1,gz1,ell1,ell_pa1,sh_str1,sh_pa1,idx1=np.loadtxt('intmpar1.txt',usecols=(0,1,2,5,6,7,8,9,10,11,12,13,14,15),unpack=True);

    if(gra1.ndim==0):
	gfld1,indx1=np.loadtxt('intmpar1.txt',dtype={'names':('gfld1','indx1'),'formats':('S13','S3')},usecols=(3,4),unpack=True);
    else:
	gfld1,indxno1=np.loadtxt('intmpar1.txt',dtype={'names':('gfld1','indxno1'),'formats':('S13','S3')},usecols=(3,4),unpack=True);
     
    zsc_q,qu,qg,qr,qi,qz=np.loadtxt(bkgqsocatalog,usecols=(3,5,6,7,8,9),unpack=True);

    idxt1,sigm1,rein1,smii1,zs1=np.loadtxt("all_qso.txt",unpack=True);
    
    rein1=rein1/pixsc;

    flxrefq=10.**(-0.4*(mi_lim_cfht_q-zpt));
    flxlimbrt=10.**(-0.4*(mi_totbrt_q-zpt));
    
    ## Set the flagg to 0, if you don't want the gravlens input files in mckq to
    ## recreated (see function srcposrnq in gensrcpos2.py )
    flagg=1;
    
    ## Output catalogs
    fp1=open('fpar1_%d.txt'%(num),'w');

    np.random.seed(num*10+29824);

    ## For bkg qso properties
    idxtt1,uid1= np.unique(idxt1, return_index=True)
    
    reinst1=np.arange(idxtt1.size)*0.;
    smi_n1=np.arange(idxtt1.size)*0.;
    zs_n1=np.arange(idxtt1.size)*0.;
    srcx1=np.arange(idxtt1.size)*0.;
    srcy1=np.arange(idxtt1.size)*0.;
    smag1=np.arange(idxtt1.size)*0.;
    imno1=np.arange(idxtt1.size)*0;
    
    smu1=np.arange(gid1.size)*0.;
    smg1=np.arange(gid1.size)*0.;
    smr1=np.arange(gid1.size)*0.;
    smi1=np.arange(gid1.size)*0.;
    smz1=np.arange(gid1.size)*0.;
    
    
    iimin=np.int(zd1.size*num/nproc);
    iimax=np.int(zd1.size*(num+1)/nproc);
    if(num==nproc-1):
	iimax=zd1.size;
    
    ## For the following loop to work, the all_xxx.txt and intmpar?.txt files
    ## need to be sorted in ascending order
    
    ## Set the range within which to match magnitudes and redshift of real quasars in
    ## order to extract colors 
   #magdiff=0.1;
   #zdiffq=0.1;
    
    ## This loop is run for each lens in the lens catalog
    for ii in range(iimin,iimax):
        ############
    ## PART 1- Select one bkg source from multiple sources 
        ############
        kkl=uid1[ii];
        if(ii==idxtt1.size-1):
            kkh=idxt1.size;
        else:
            kkh=uid1[ii+1];
    
        ## For each source, extract source position, flux of the 2nd brightest image, total
        ## magnification of the lensed source and number of lensed images
        srcxt,srcyt,smagt,sumt,imtno=scp.srcposrnq(rein1[kkl:kkh],ell1[ii],ell_pa1[ii],sh_str1[ii],sh_pa1[ii],ii+1,flagg,int(gid1[ii]),num);
        
        cnt1=np.arange(kkh-kkl)*0;
        jj=0;
        kk=kkl;
        for kk in range(kkl,kkh):
            hh=kkl-kk;
            flx2=10.**(-0.4*(smii1[kk]-zpt)) * smagt[hh];
            flxall=10.**(-0.4*(smii1[kk]-zpt)) * sumt[hh];
            ## Accept each source, if the 2nd brightest lensed image and sum of flux of all lensed
            ## images is above the set limits
            if(flx2>flxrefq and flxall<flxlimbrt):
                cnt1[jj]=kk;
                jj=jj+1;
        
        ## Choose one source randomly from the sources which satisfy the flux
        ## limits
        if(jj>0):
            qq=np.random.randint(0,jj);
        else:
            print "No source with 2nd brightest image above mlim=",mi_lim_cfht_q,"for lens id:",int(gid1[ii]);
            continue;
        nq=cnt1[qq]-kkl;
        smi_n1[ii]=smii1[cnt1[qq]];
        zs_n1[ii]=zs1[cnt1[qq]];
        reinst1[ii]=rein1[cnt1[qq]];
        srcx1[ii]=srcxt[nq]
        srcy1[ii]=srcyt[nq]
        smag1[ii]=smagt[nq]
        imno1[ii]=imtno[nq]
        sys.stdout.flush();
        
        ############
        ## PART 2- Extract qso colors from real qso catalog
            ############
        ll=0;
        indx_n=np.arange(qu.size)*0;
        for jj in range(qu.size):
            if(abs(smi_n1[ii]-qi[jj])<=magdiff and abs(zs_n1[ii]-zsc_q[jj])<=zdiffq_strict):
                indx_n[ll]=jj;
                ll=ll+1;
        if(ll==0):
            for jj in range(qu.size):
                if (abs(zs_n1[ii]-zsc_q[jj])<=zdiffq_relaxed and qg[jj]-qi[jj]<gicolorq_relaxed):
                    indx_n[ll]=jj;
                    ll=ll+1;

        indx_n=indx_n[0:ll];
        nn=np.random.randint(0,ll);
        ncnt1=indx_n[nn];
        
        ## Use qso colors only (not the magnitudes) from the catalog
        smi1[ii]=smi_n1[ii];
        smu1[ii]=smi_n1[ii]+qu[ncnt1]-qi[ncnt1];
        smg1[ii]=smi_n1[ii]+qg[ncnt1]-qi[ncnt1];
        smr1[ii]=smi_n1[ii]+qr[ncnt1]-qi[ncnt1];
        smz1[ii]=smi_n1[ii]+qz[ncnt1]-qi[ncnt1];
        
	## Save catalog with all the lens+source parameters
	fp1.write('%d %f %f %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d \n'%(gid1[ii],gra1[ii],gdec1[ii],gfld1[ii],indxno1[ii],zd1[ii],gu1[ii],gg1[ii],gr1[ii],gi1[ii],gz1[ii],ell1[ii],ell_pa1[ii],reinst1[ii],sh_str1[ii],sh_pa1[ii],srcx1[ii],srcy1[ii],smu1[ii],smg1[ii],smr1[ii],smi1[ii],smz1[ii],zs_n1[ii],smag1[ii],imno1[ii]));
    
	## Generate qso lensed images
	gn.genimg_gq(reinst1[ii],ell1[ii],ell_pa1[ii],sh_str1[ii],sh_pa1[ii],smu1[ii],smg1[ii],smr1[ii],smi1[ii],smz1[ii],srcx1[ii],srcy1[ii],zpt,gid1[ii],indxno1[ii],num);
    
    fp1.close();    

## Run this code faster by specifying Nproc (no. of processors)
jobs=[];
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc));
    jobs.append(p);
    p.start();

print  "#####################################################################";
