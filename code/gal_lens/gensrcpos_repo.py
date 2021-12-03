#!/usr/bin/env python
from math import *
from subprocess import call
import numpy as np
import sys
import random
import os    
from input_qg import *

np.random.seed(72398);

## Extract source position and some lensed image properties for quasars
def srcposrnq(reinst,ell,ell_pa,sh_str,sh_pa,glno,flag,gid,num):
    limr=reinst.size;
    if(flag):
	fp=open('mckq/gllens%d.in'%(glno),'w');

	np.savetxt(fp,['set rscale=15'],fmt='%s');
	np.savetxt(fp,['set ngrid1=40'],fmt='%s');
	np.savetxt(fp,['set ngrid2=40'],fmt='%s');
	np.savetxt(fp,['set verbose=0'],fmt='%s');
	np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');

	for ii in range(limr):
	    np.savetxt(fp,['startup 1 1'],fmt='%s');
	    fp.write('alpha %f 0 0 %f %f %f %f 0 0 1 \n' %(reinst[ii],ell,ell_pa,sh_str,sh_pa));
	    np.savetxt(fp,[' 0 0 0 0 0 0 0 0 0 0'],fmt='%s');
	    fp.write("mock1 mckq/out%d_%d 2000 2 \n \n" %(glno,ii+1));

	np.savetxt(fp,['quit'],fmt='%s');
	fp.close();
     
    print "Running mckq/gllens%d.in on %d"%(glno,num);
    call('%s mckq/gllens%d.in > mckq/gltmpout_%d_%d'%(lenscode,glno,glno,num),shell=1);
    
    ii=0;
    srcx=np.arange(limr)*0.;
    srcy=np.arange(limr)*0.;
    smag=np.arange(limr)*0.;
    summ=np.arange(limr)*0.;
    imno=np.arange(limr)*0;
     
    for ii in range(limr):
        ip=ii+1;
	imgnos=np.random.uniform(0,1);
	## Loop to select a source that has 2 or more lensed images such that a
	## double is chosen 50% of the times
	cond=1;
	while(cond<3):
	    if(imgnos<0.5):
		call("awk -f scr2img mckq/out'%d'_'%d' > mckq/tp_%d"%(glno,ip,glno),shell=True);
		## if no doubly images sources (unlikely to happen), try sources with no. of lensed
		## images >2
		if( os.path.getsize("mckq/tp_%d"%(glno))<=0):
		    imgnos=0.6;
		    cond=cond+1;
		else:
		    rwno=np.loadtxt("mckq/tp_%d"%(glno),unpack=True);
		    cond=10;
	    else:    
		call("awk \'{if($NF~/images/ && $1>2) print NR}\' mckq/out'%d'_'%d' > mckq/tp_%d"%(glno,ip,glno),shell=True);
		## if no sources with no. of lensed images >2, try to find
		## doubly imaged sources
		if( os.path.getsize("mckq/tp_%d"%(glno))<=0):
		    imgnos=0.4;
		    cond=cond+1;
		else:
		    rwno=np.loadtxt("mckq/tp_%d"%(glno),unpack=True);
		    cond=10;
        if(cond==3):
	    print "No lensed images found for lens ",ii,"with id",gid,"..continuing..";
	    continue;
        else: 
	    ## More than one eligible sources found 
	    if(rwno.ndim>0):    
		indx=np.random.randint(0,rwno.size);
		rwn=int(rwno[indx]);
		
		filename=('mckq/out%d_%d'%(glno,ip));
	        
		## Extract source position and no. of lensed images
		call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d \\n\",$1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s > mckq/tmp_%d"%(filename,glno),shell=1);
		## Extract source position, no. of lensed images, row no. of the
		## first lensed image and lens galaxy id
		call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d %d '+"%d"%(gid)+' \\n\",$2,NR+1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s >> mckq/srcinfo_%d"%(filename,glno),shell=1);
		xx,yy,jj=np.loadtxt('mckq/tmp_%d'%(glno),unpack=True);
		call("awk \'{if(NR>'%d' && NR<='%d') print $3}\' '%s'  >  mckq/tmp2_%d" % (rwn,rwn+jj,filename,glno),shell=True);
		magn=np.loadtxt('mckq/tmp2_%d'%(glno),unpack=True);
		magnsrt=sorted(np.abs(magn),reverse=1);
		magsum=sum(np.abs(magn));
		
		srcx[ii]=xx;
		srcy[ii]=yy;
		smag[ii]=magnsrt[1]; 
		summ[ii]=magsum; 
		imno[ii]=int(jj);
		
	    ## One eligible source found 
	    else:
		rwn=int(rwno);
		
		filename=('mckq/out%d_%d'%(glno,ip));
		call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d \\n\",$1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s > mckq/tmp_%d"%(filename,glno),shell=1);
		call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d %d '+"%d"%(gid)+' \\n\",$2,NR+1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s >> mckq/srcinfo_%d"%(filename,glno),shell=1);
		xx,yy,jj=np.loadtxt('mckq/tmp_%d'%(glno),unpack=True);
		srcx[ii]=xx;
		srcy[ii]=yy;
		imno[ii]=int(jj);

		call("awk \'{if(NR>'%d' && NR<='%d') print $3}\' mckq/out'%d'_'%d' >  mckq/tmp2_%d" % (rwn,rwn+jj,glno,ip,glno),shell=True);
		magn=np.loadtxt('mckq/tmp2_%d'%(glno),unpack=True);
		magnsrt=sorted(np.abs(magn),reverse=1);
		magsum=sum(np.abs(magn));
		smag[ii]=magnsrt[1]; 
		summ[ii]=magsum; 

    return srcx,srcy,smag,summ,imno;


## Extract source position and some lensed image properties for galaxies
def srcposrng(reinst,ell,ell_pa,sh_str,sh_pa,glno,flag,gid,num):
    limr=reinst.size;
    if(flag):
	fp=open('mckg/gllens%d.in'%(glno),'w');

	np.savetxt(fp,['set rscale=15'],fmt='%s');
	np.savetxt(fp,['set ngrid1=40'],fmt='%s');
	np.savetxt(fp,['set ngrid2=40'],fmt='%s');
	np.savetxt(fp,['set verbose=0'],fmt='%s');
	np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');

	for ii in range(limr):
	    np.savetxt(fp,['startup 1 1'],fmt='%s');
	    fp.write('alpha %f 0 0 %f %f %f %f 0 0 1 \n' %(reinst[ii],ell,ell_pa,sh_str,sh_pa));
	    np.savetxt(fp,[' 0 0 0 0 0 0 0 0 0 0'],fmt='%s');
	    fp.write("mock1 mckg/out%d_%d 2000 2 \n \n" %(glno,ii+1));

	np.savetxt(fp,['quit'],fmt='%s');
	fp.close();
     
	print "Running mckg/gllens%d.in on %d"%(glno,num);
	call('%s mckg/gllens%d.in > mckg/gltmpout_%d_%d'%(lenscode,glno,glno,num),shell=1);
    
    ii=0;
    srcx=np.arange(limr)*0.;
    srcy=np.arange(limr)*0.;
    smag=np.arange(limr)*0.;
    summ=np.arange(limr)*0.;
    imno=np.arange(limr)*0;
     
    for ii in range(limr):
        ip=ii+1;
	imgnos=np.random.uniform(0,1);
	## Loop to select a source that has 2 or more lensed images such that a
	## double is chosen 50% of the times
	cond=1;
	while(cond<3):
	    if(imgnos<0.5):
		call("awk -f scr2img mckg/out'%d'_'%d' > mckg/tp_%d"%(glno,ip,glno),shell=True);
		## if no doubly images sources (unlikely to happen), try sources with no. of lensed
		## images >2
		if( os.path.getsize("mckg/tp_%d"%(glno))<=0):
		    imgnos=0.6;
		    cond=cond+1;
		else:
		    rwno=np.loadtxt("mckg/tp_%d"%(glno),unpack=True);
		    cnt=0;
		    cond=10;
	    else:    
		call("awk \'{if($NF~/images/ && $1>2) print NR}\' mckg/out'%d'_'%d' > mckg/tp_%d"%(glno,ip,glno),shell=True);
		## if no sources with no. of lensed images >2, try to find
		## doubly imaged sources
		if( os.path.getsize("mckg/tp_%d"%(glno))<=0):
		    imgnos=0.4;
		    cond=cond+1;
		else:
		    rwno=np.loadtxt("mckg/tp_%d"%(glno),unpack=True);
		    cnt=0;
		    cond=10;
        if(cond==3):
	    print "No lensed images found for lens ",ii,"with id",gid,"..continuing..";
	    continue;
        
	restart_c=True;
        ## Repeat this loop, if condition on excluding ER configuration fail for the current iteration
	while restart_c:
	    ## Avoid selecting a large fraction of Einstein Ring configurations, by excluding highly aligned sources with the lenses
	    condition=True;
	    while condition:
		## if jj>2 fails after 5 attempts
		if(cnt>5):
		    call("awk -f scr2img mckg/out'%d'_'%d' > mckg/tp_%d"%(glno,ip,glno),shell=True);
	
		    if( os.path.getsize("mckg/tp_%d"%(glno))<=0):
			break;
		    else:
			rwno=np.loadtxt("mckg/tp_%d"%(glno),unpack=True);
			cnt=0;
		
		## More than one eligible sources found 
		if(rwno.ndim>0):    
		    indx=np.random.randint(0,rwno.size);
		    rwn=int(rwno[indx]);
		    
		    filename=('mckg/out%d_%d'%(glno,ip));
		    print "gllens%d.in"%(glno),filename,rwn;
		    
		    ## Extract source position and no. of lensed images
		    call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d \\n\",$1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s > mckg/tmp_%d"%(filename,glno),shell=1);
		    ## Extract source position, no. of lensed images, row no. of the
		    ## first lensed image and lens galaxy id
		    call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d %d '+"%d"%(gid)+' \\n\",$2,NR+1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s >> mckg/srcinfo_%d"%(filename,glno),shell=1);
		    xx,yy,jj=np.loadtxt('mckg/tmp_%d'%(glno),unpack=True);
		    
		    if( (abs(xx)>=0.1 and abs(yy)>=0.1)  or (cnt>5 and jj==2)):
			call("awk \'{if(NR>'%d' && NR<='%d') print $3}\' '%s'  >  mckg/tmp2_%d" % (rwn,rwn+jj,filename,glno),shell=True);
			magn=np.loadtxt('mckg/tmp2_%d'%(glno),unpack=True);
			magnsrt=sorted(np.abs(magn),reverse=1);
			magsum=sum(np.abs(magn));
			
			srcx[ii]=xx;
			srcy[ii]=yy;
			smag[ii]=magnsrt[1]; 
			summ[ii]=magsum; 
			imno[ii]=int(jj);
			condition=False;
			restart_c=False;
		    else:
			condition=True; 
			cnt=cnt+1;

		    if(cnt>5 and jj>2):
                        ## Try extracting xx,yy at least 5 times before accepting an ER-like configuration
			jj=jj-1;
			restart_c=True;
			condition=False;
			
		else:
		    rwn=int(rwno);
		    
		    filename=('mckg/out%d_%d'%(glno,ip));
		    print "gllens%d.in"%(glno),filename,rwn;
		    
		    ## Extract source position and no. of lensed images
		    call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d \\n\",$1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s > mckg/tmp_%d"%(filename,glno),shell=1);
		    ## Extract source position, no. of lensed images, row no. of the
		    ## first lensed image and lens galaxy id
		    call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d %d '+"%d"%(gid)+' \\n\",$2,NR+1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s >> mckg/srcinfo_%d"%(filename,glno),shell=1);
		    xx,yy,jj=np.loadtxt('mckg/tmp_%d'%(glno),unpack=True);
		    srcx[ii]=xx;
		    srcy[ii]=yy;
		    imno[ii]=int(jj);

		    call("awk \'{if(NR>'%d' && NR<='%d') print $3}\' mckg/out'%d'_'%d' >  mckg/tmp2_%d" % (rwn,rwn+jj,glno,ip,glno),shell=True);
		    magn=np.loadtxt('mckg/tmp2_%d'%(glno),unpack=True);
		    magnsrt=sorted(np.abs(magn),reverse=1);
		    magsum=sum(np.abs(magn));
		    smag[ii]=magnsrt[1]; 
		    summ[ii]=magsum; 
		    
		    condition=False;
		    restart_c=False;
    return srcx,srcy,smag,summ,imno;


## Calculate size of the source assuming size-luminosity relation
## using Bernardi et al. 2003 Eqns given in Oguri 2006
def srcsize(mapp,zsrc,pixsc):
    ## mapp is in g-band
    Dlum=cc.Dlofz(zsrc)/p.hval;
    Mabs =mapp-5*log10(Dlum)-25;
    Lum_src=10**(-0.4*(Mabs-5.48));
    Da=cc.Daofz(zsrc)/p.hval*1.e3;

    Lrat= Lum_src/10**10.2;
    Reff= 10**0.52*Lrat**(2./3.) * 1./(1+zsrc)**2;
    ## Half-light radius converting from kpc->radian->arcsec->pix
    return (Reff/Da)*(180.*3600/pi/pixsc);
