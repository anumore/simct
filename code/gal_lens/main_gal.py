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
## mkcatalog.py and Background galaxy catalog which has Redshift and Magnitude 
##
## PURPOSE:
## Use the lens properties to choose a single source per lens, calculate 
## source position based on the limits on fluxes of the lensed images, assign
## colors and generate input files which will produce lensed galaxy images 
###########################################################

## Galaxies as background sources

def worker(num,nproc):
    ## Read the lens galaxy catalog, real bkg galaxy color catalog and 
    ## bkg galaxy catalog generated for each potential lens 
    gid0,gra0,gdec0,zd0,gu0,gg0,gr0,gi0,gz0,ell0,ell_pa0,sh_str0,sh_pa0,idx0=np.loadtxt('intmpar0.txt',usecols=(0,1,2,5,6,7,8,9,10,11,12,13,14,15),unpack=True);

    if(gra0.ndim==0):
        gfld0,indx0=np.loadtxt('intmpar0.txt',dtype={'names':('gfld0','indx0'),'formats':('S13','S3')},usecols=(3,4),unpack=True);
    else:
        gfld0,indxno0=np.loadtxt('intmpar0.txt',dtype={'names':('gfld0','indxno0'),'formats':('S13','S3')},usecols=(3,4),unpack=True);
     
    zsc_g,zsc_mn,zsc_mx,su,sg,sr,si,sz=np.loadtxt(bkggalcatalog,usecols=(5,6,7,8,9,10,11,12),unpack=True);

    idxt0,sigm0,rein0,smii0,zs0=np.loadtxt("all_gal_mod.txt",unpack=True);
    
    rein0=rein0/pixsc;

    flxrefg=10.**(-0.4*(mi_lim_cfht_g-zpt));
    flxlimbrt=10.**(-0.4*(mi_totbrt_g-zpt));
    
    ## Set the flagg to 0, if you don't want the gravlens input files in mckg to
    ## recreated (see function srcposrnq in gensrcpos2.py )
    flagg=1;
    
    ## Output catalogs
    fp0=open('fpar0_%d.txt'%(num),'w');

    np.random.seed(num*10+29824);

    ## For bkg gal properties
    idxtt0,uid0= np.unique(idxt0, return_index=True);

    reinst0=np.arange(idxtt0.size)*0.;
    smi_n0=np.arange(idxtt0.size)*0.;
    zs_n0=np.arange(idxtt0.size)*0.;
    srcx0=np.arange(idxtt0.size)*0.;
    srcy0=np.arange(idxtt0.size)*0.;
    smag0=np.arange(idxtt0.size)*0.;
    imno0=np.arange(idxtt0.size)*0;
    
    smu0=np.arange(gid0.size)*0.;
    smg0=np.arange(gid0.size)*0.;
    smr0=np.arange(gid0.size)*0.;
    smi0=np.arange(gid0.size)*0.;
    smz0=np.arange(gid0.size)*0.;
    
    ell_s=np.arange(gid0.size)*0.;
    ell_pa_s=np.arange(gid0.size)*0.;
    reff_s=np.arange(gid0.size)*0.;
    
    iimin=np.int(zd0.size*num/nproc);
    iimax=np.int(zd0.size*(num+1)/nproc);
    if(num==nproc-1):
        iimax=zd0.size;
    
    ## For the following loop to work, the all_xxx.txt and intmpar?.txt files
    ## need to be sorted in ascending order
    
    ## Set the range within which to match magnitudes and redshift of real galaxies 
    ## in order to extract colors 
    #magdiff=0.1;
    zdiffg=zsc_mx-zsc_mn;
    
    ## This loop is run for each lens in the lens catalog
    for ii in range(iimin,iimax):
        ############
	## PART 1- Select one bkg source from multiple sources 
        ############
        kkl=uid0[ii];
        if(ii==idxtt0.size-1):
            kkh=idxt0.size;
        else:
            kkh=uid0[ii+1];
	
	## For each source, extract source position, flux of the 2nd brightest image, total
	## magnification of the lensed source and number of lensed images
        srcxt,srcyt,smagt,sumt,imtno=scp.srcposrng(rein0[kkl:kkh],ell0[ii],ell_pa0[ii],sh_str0[ii],sh_pa0[ii],ii+1,flagg,int(gid0[ii]),num);
    
        cnt0=np.arange(kkh-kkl)*0;
        jj=0;
        kk=kkl;
        for kk in range(kkl,kkh):
            hh=kk-kkl;
            flx2=10.**(-0.4*(smii0[kk]-zpt)) * smagt[hh];
            flxall=10.**(-0.4*(smii0[kk]-zpt)) * sumt[hh];
        ## Accept each source, if the 2nd brightest lensed image and sum of flux of all lensed
        ## images is above the set limits
            if(flx2>flxrefg and flxall<flxlimbrt):
                cnt0[jj]=kk;
                jj=jj+1;
        
        ## Choose one source randomly from the sources which satisfy the flux
        ## limits
        if(jj>0):
            qq=np.random.randint(0,jj);
        else:
            print "No source with 2nd brightest image above mlim=",mi_lim_cfht_g,"for lens id:",int(gid0[ii]);
            continue;
        nq=cnt0[qq]-kkl;
        smi_n0[ii]=smii0[cnt0[qq]];
        zs_n0[ii]=zs0[cnt0[qq]];
        reinst0[ii]=rein0[cnt0[qq]];
        srcx0[ii]=srcxt[nq];
        srcy0[ii]=srcyt[nq];
        smag0[ii]=smagt[nq];
        imno0[ii]=imtno[nq];
    
        ############
	## PART 2- Extract gal colors from real gal catalog
        ############
        ll=0;
        indx_n=np.arange(su.size)*0;
        for jj in range(su.size):
            if (abs(smi_n0[ii]-si[jj])<=magdiff and abs(zs_n0[ii]-zsc_g[jj])<=zdiffg_strict*zdiffg[jj] and sg[jj]-si[jj]<gicolorg_strict):
                indx_n[ll]=jj;
                ll=ll+1;
        if(ll==0):
            for jj in range(su.size):
                if (abs(zs_n0[ii]-zsc_g[jj])<=zdiffg_relaxed*zdiffg[jj] and sg[jj]-si[jj]<gicolorg_relaxed):
                    indx_n[ll]=jj;
                    ll=ll+1;

        indx_n=indx_n[0:ll];
        nn=np.random.randint(0,ll);
        ncnt0=indx_n[nn];
        ## Use gal colors only (not the magnitudes) from the catalog
        smi0[ii]=smi_n0[ii];
        smu0[ii]=smi_n0[ii]+su[ncnt0]-si[ncnt0];
        smg0[ii]=smi_n0[ii]+sg[ncnt0]-si[ncnt0];
        smr0[ii]=smi_n0[ii]+sr[ncnt0]-si[ncnt0];
        smz0[ii]=smi_n0[ii]+sz[ncnt0]-si[ncnt0];

        ## Generate random ellipticities, PA and size for the background source
        ell_s[ii]=np.random.uniform(bkggal_ell_low,bkggal_ell_high);
        ell_pa_s[ii]=np.random.uniform(bkggal_ellpa_low,bkggal_ellpa_high);
        reff_s[ii]=scp.srcsize(smg0[ii],zs_n0[ii],pixsc);

	## Save catalog with all the lens+source parameters
        fp0.write('%d %f %f %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %f\n'%(gid0[ii],gra0[ii],gdec0[ii],gfld0[ii],indxno0[ii],zd0[ii],gu0[ii],gg0[ii],gr0[ii],gi0[ii],gz0[ii],ell0[ii],ell_pa0[ii],reinst0[ii],sh_str0[ii],sh_pa0[ii],srcx0[ii],srcy0[ii],smu0[ii],smg0[ii],smr0[ii],smi0[ii],smz0[ii],zs_n0[ii],smag0[ii],imno0[ii],ell_s[ii],ell_pa_s[ii],reff_s[ii]));

        ## Generate lensmodel inp files
        gn.genimg_gg(reinst0[ii],ell0[ii],ell_pa0[ii],sh_str0[ii],sh_pa0[ii],smu0[ii],smg0[ii],smr0[ii],smi0[ii],smz0[ii],srcx0[ii],srcy0[ii],ell_s[ii],ell_pa_s[ii],reff_s[ii],zpt,gid0[ii],indxno0[ii],num);
	
    fp0.close();    

## Run this code faster by specifying Nproc (no. of processors)
jobs=[];
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc));
    jobs.append(p);
    p.start();

print  "#####################################################################";
