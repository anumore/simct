#!/usr/bin/env python
from math import *
import numpy as np
import sys
import genimg as gn
import srclensprop as slp
from subprocess import call
import multiprocessing
import time
from input_g import *

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
    ## and lens catalog with info on BCG+members
    gid0,gra0,gdec0,zd0,gu0,gg0,gr0,gi0,gz0,ell0,ell_pa0=np.loadtxt('intmpar.txt',usecols=(0,1,2,5,6,7,8,9,10,11,12),unpack=True);

    if(gra0.ndim==0):
	gfld0,indx0=np.loadtxt('intmpar.txt',dtype={'names':('gfld0','indx0'),'formats':('S13','S3')},usecols=(3,4),unpack=True);
    else:
	gfld0,indxno0=np.loadtxt('intmpar.txt',dtype={'names':('gfld0','indxno0'),'formats':('S13','S3')},usecols=(3,4),unpack=True);
     
    zsc_g,zsc_mn,zsc_mx,su,sg,sr,si,sz=np.loadtxt(bkggalcatalog,usecols=(5,6,7,8,9,10,11,12),unpack=True);

    nwmid,idxt0,ouid,nsrc,smii0,zs0=np.loadtxt("all_gal.txt",usecols=(0,1,2,4,5,6),unpack=True);
    
    mid,gid,kappa1,rscale,gposx,gposy,zds,sig,ell,ell_pa=np.loadtxt("lenscat.txt",unpack=True);

    flxrefg=10.**(-0.4*(mi_lim_cfht_g-zpt));
    flxlimbrt=10.**(-0.4*(mi_totbrt_g-zpt));
    
    ## Set the flagg to 0, if you don't want the gravlens input files in mckg to
    ## recreated (see function srcposrnq in gensrcpos2.py )
    flagg=1;
    
    ## Calculate Einstein radius for each BCG + member and kappa for every lens
    nwid,nouid=np.unique(ouid,return_index=True);
    kappa=np.arange(zs0.size)*0.;
    reinst=np.arange(mid.size*zs0.size).reshape(zs0.size,mid.size)*0.;
    jj=0;
    previd=-99
    for kk in range(zs0.size):
	if(ouid[kk]!=previd):
	    xxl=int(nwid[jj]);
	    if(jj==nwid.size-1):
	        xxh=zds.size;
	    else:
	        xxh=int(nwid[jj+1]);
	    jj=jj+1;   
	reint,drat=slp.getreinst2(zds[xxl],zs0[kk],sig[xxl:xxh])
	for xx in range(reint.size):
	    reinst[kk][xx]=reint[xx]*(206264.8/pixsc);
	kappa[kk]=kappa1[xxl]/drat;
	previd=ouid[kk];
    
    ## For bkg gal properties
    smi_n0=np.arange(nwid.size)*0.;
    zs_n0=np.arange(nwid.size)*0.;
    srcx0=np.arange(nwid.size)*0.;
    srcy0=np.arange(nwid.size)*0.;
    smag0=np.arange(nwid.size)*0.;
    imno0=np.arange(nwid.size)*0;
    kappa0=np.arange(nwid.size)*0.;
    reinst0=np.arange(mid.size*nwid.size).reshape(nwid.size,mid.size);
     
    smu0=np.arange(gid0.size)*0.;
    smg0=np.arange(gid0.size)*0.;
    smr0=np.arange(gid0.size)*0.;
    smi0=np.arange(gid0.size)*0.;
    smz0=np.arange(gid0.size)*0.;

    ell_s=np.arange(gid0.size)*0.;
    ell_pa_s=np.arange(gid0.size)*0.;
    reff_s=np.arange(gid0.size)*0.;
    
    print "Extracting bkg gal properties";
    sys.stdout.flush();
    

    iimin=np.int(zd0.size*num/nproc);
    iimax=np.int(zd0.size*(num+1)/nproc);
    if(num==nproc-1):
	iimax=zd0.size;

    ## Output catalogs
    fp0=open('fpar_%d.txt'%(num),'w');
    
    np.random.seed(num*10+29924);
    
    ## For the following loop to work, the all_xxx.txt and intmpar?.txt files
    ## need to be sorted in ascending order
    
    ## Set the range within which to match magnitudes and redshift of real galaxies 
    ## in order to extract colors 
    ##magdiff=0.1;
    zdiffg=zsc_mx-zsc_mn;
    
    for ii in range(iimin,iimax):
        ############
	## PART 1- Select one bkg source from multiple sources 
        ############
	## kk is for sources per lens or just every lens
	## xx is for group members per lens
	kk=int(nouid[ii]);	
	if(ouid[kk]!=ouid[int(nouid[ii-1])] or kk==0):
	    kkl=kk;
	    kkh=kk+int(nsrc[kk]);
	    xxl=int(nwid[ii]);
	    if(ii==nwid.size-1):
		xxh=zds.size;
	    else:
		xxh=int(nwid[ii+1]);

	## For each source, extract source position, flux of the 2nd brightest image, total
	## magnification of the lensed source and number of lensed images
	srcxt,srcyt,smagt,sumt,imtno=gn.srcpos(kappa[kkl:kkh],rscale[xxl],reinst[kkl:kkh],gposx[xxl:xxh],gposy[xxl:xxh],ell[xxl:xxh],ell_pa[xxl:xxh],ii+1,int(nsrc[kkl]),int(gid[xxl]),num);
       
	cnt=np.arange(kkh-kkl)*0;
        jj=0;
        kk=kkl;
        for kk in range(kkl,kkh):
            hh=kk-kkl;
            sys.stdout.flush();
            flx2=10.**(-0.4*(smii0[kk]-zpt)) * smagt[hh];
            flxall=10.**(-0.4*(smii0[kk]-zpt)) * sumt[hh];
	    ## Accept each source, if the 2nd brightest lensed image and sum of flux of all lensed
	    ## images is above the set limits
            if(flx2>flxrefg and flxall<flxlimbrt):
        	cnt[jj]=kk;
        	jj=jj+1;
	
	## Choose one source randomly from the sources which satisfy the flux
	## limits
	if(jj>0):
	    qq=np.random.randint(0,jj);
	else:
	    print "No source with 2nd brightest image above mlim=",mi_lim_cfht_g,"for lens id:",int(gid0[ii]);
	    continue;
        nq=cnt[qq]-kkl;
        smi_n0[ii]=smii0[cnt[qq]];
        zs_n0[ii]=zs0[cnt[qq]];
        kappa0[ii]=kappa[cnt[qq]];
        reinst0[ii]=reinst[cnt[qq]];
        srcx0[ii]=srcxt[nq]
        srcy0[ii]=srcyt[nq]
        smag0[ii]=smagt[nq]
        imno0[ii]=imtno[nq]

        ############
	## PART 2- Extract gal colors from real gal catalog
        ############
	ll=0;
        indx_n=np.arange(su.size)*0;
        for jj in range(su.size):
            if (abs(smi_n0[ii]-si[jj])<=magdiff and abs(zs_n0[ii]-zsc_g[jj])<=zdiffg_strict*zdiffg[jj] and sg[jj]-si[jj]<gicolorg_strict ):
		indx_n[ll]=jj;
        	ll=ll+1;
        if(ll==0):
            for jj in range(su.size):
        	if (abs(zs_n0[ii]-zsc_g[jj])<=zdiffg_relaxed*zdiffg[jj] and sg[jj]-si[jj]<gicolorg_relaxed):
        	    indx_n[ll]=jj;
        	    ll=ll+1;
        
        indx_n=indx_n[0:ll];
        nn=np.random.randint(0,ll);
        ncnt=indx_n[nn];
        ## Use gal colors only (not the magnitudes) from the catalog
        smi0[ii]=smi_n0[ii];
        smu0[ii]=smi_n0[ii]+su[ncnt]-si[ncnt];
        smg0[ii]=smi_n0[ii]+sg[ncnt]-si[ncnt];
        smr0[ii]=smi_n0[ii]+sr[ncnt]-si[ncnt];
        smz0[ii]=smi_n0[ii]+sz[ncnt]-si[ncnt];

        ## Generate random ellipticities, PA and size for the background source
        ell_s[ii]=np.random.uniform(bkggal_ell_low,bkggal_ell_high);
        ell_pa_s[ii]=np.random.uniform(bkggal_ellpa_low,bkggal_ellpa_high);
        reff_s[ii]=gn.srcsize(smg0[ii],zs_n0[ii],pixsc);

	## Save catalog with all the lens+source parameters
        fp0.write('%d %f %f %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %f %f %f\n'%(gid0[ii],gra0[ii],gdec0[ii],gfld0[ii],indxno0[ii],zd0[ii],gu0[ii],gg0[ii],gr0[ii],gi0[ii],gz0[ii],ell0[ii],ell_pa0[ii],srcx0[ii],srcy0[ii],smu0[ii],smg0[ii],smr0[ii],smi0[ii],smz0[ii],zs_n0[ii],smag0[ii],imno0[ii],ell_s[ii],ell_pa_s[ii],reff_s[ii]));

        print "Generating bkg gal lensmodel inp files"
        gn.genimg_gg(kappa0[ii],rscale[xxl],reinst0[ii],gposx[xxl:xxh],gposy[xxl:xxh],ell[xxl:xxh],ell_pa[xxl:xxh],smu0[ii],smg0[ii],smr0[ii],smi0[ii],smz0[ii],srcx0[ii],srcy0[ii],ell_s[ii],ell_pa_s[ii],reff_s[ii],zpt,gid0[ii],indxno0[ii],int(nsrc[kkl]),num);
    
    print "Saving all the lens+source properties"
    sys.stdout.flush();
    fp0.close();    

## Run this code faster by specifying Nproc (no. of processors)
jobs=[];
for i in range(Nproc):
    p = multiprocessing.Process(target=worker,args=(i,Nproc));
    jobs.append(p);
    p.start();

print  "#####################################################################";
