#!/usr/bin/env python
from math import *
import numpy as np
from subprocess import call
import os
import pyfits
from scipy.interpolate import interp1d
import srclensprop as slp
from input_g import *


## Spatial scales are in pixels and flux is converted to counts/sec

init_cs_spl=0;
cross_sect_spl=0;

## Calculate lens cross-section
def getcsarea(kappa,rscale,gposx,gposy,reinst,ell,ell_pa,ii,pixsc,num):
    ## Initialize lcno, get kappa and gposx,gposy
    ## create the gravlens file for each lens
    call('rm -f ./mck1/cc%d\n'%(ii),shell=1);
    lcno=reinst.size+1;
    fvar='';
    fvar += ('mck1/gclens%s_%d.in' % (ii,num));
    fp=open(fvar,'w');

    np.savetxt(fp,['set rscale=15'],fmt='%s');
    np.savetxt(fp,['set ngrid1=40'],fmt='%s');
    np.savetxt(fp,['set ngrid2=40'],fmt='%s');
    np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');

    fp.write('startup %d 1 \n'%(lcno));
    fp.write('nfw %f 0 0 0 0 0 0 %f 0 0 \n' %(kappa,rscale));
    for mm in range(reinst.size):
	fp.write('alpha %f %f %f %f %f 0 0 0 0 1 \n' %(reinst[mm],gposx[mm],gposy[mm],ell[mm],ell_pa[mm]));
    for mm in range(lcno): 
	np.savetxt(fp,[' 0 0 0 0 0 0 0 0 0 0 \n'],fmt='%s');

    fp.write('plotcrit ./mck1/cc%d\n'%(ii));

    np.savetxt(fp,['quit'],fmt='%s');
    fp.close();
        
    ## Calculate the area covered by caustics
    cmd_str=lenscode+" "+fvar+" >> mck1/glout  ";
    os.system(cmd_str);

    u1,v1=np.loadtxt('mck1/cc%d'%(ii),usecols=(2,3),unpack=True);
    umin=min(u1);
    umax=max(u1);
    vmin=min(v1);
    vmax=max(v1);
    csarea=(umax-umin) * (vmax-vmin);
    
    ## units of steradian
    return csarea*(pixsc/206264.8)**2;

## Initialize the spline for lens cross-section
def init_crosssect(kappa1,rscale,vdisp,gpx,gpy,ell,ellpa,zred,cnt,pixsc,num,zzmin,zzmax):
    global cross_sect_spl,init_cs_spl;
    zz=np.arange(zzmin,zzmax+(zzmax-zzmin)/10,(zzmax-zzmin)/10);
    csarea=zz*0.;
    for ii in range(zz.size):
	reinst,Drat=slp.getreinst2(zred,zz[ii],vdisp);
	kappa=kappa1/Drat;
	csarea[ii]=getcsarea(kappa,rscale,gpx,gpy,reinst*(206264.8/pixsc),ell,ellpa,cnt,pixsc,num);
    
    cross_sect_spl=interp1d(zz,csarea,kind='cubic');
    init_cs_spl=1;
    #print "Initialized cross-section spline";
    

## Calculate positions of sources and corresponding fluxes of the lensed images
## for each lens
def srcpos(kappa,rscale,reinst,gposx,gposy,ell,ell_pa,glno,nsrc,gid,num):
    ## lcno is the number of lines for the lens components, nfw+iso+iso...
    ## no. of lines to extract from an existing gravlens model file
    lcno=gposx.size+1;
    nline=(2*(lcno)+3);
    fp=open('mckg/gllens%d.in'%(glno),'w');

    np.savetxt(fp,['set rscale=15'],fmt='%s');
    np.savetxt(fp,['set ngrid1=40'],fmt='%s');
    np.savetxt(fp,['set ngrid2=40'],fmt='%s');
    np.savetxt(fp,['set verbose=0'],fmt='%s');
    np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');

    for jj in range(nsrc):
	fp.write('startup %d 1 \n'%(lcno));
	fp.write('nfw %f 0 0 0 0 0 0 %f 0 0 \n' %(kappa[jj],rscale));
	for mm in range(gposx.size):
	    fp.write('alpha %f %f %f %f %f 0 0 0 0 1 \n' %(reinst[jj][mm],gposx[mm],gposy[mm],ell[mm],ell_pa[mm]));
	for mm in range(lcno): 
	    np.savetxt(fp,['0 0 0 0 0 0 0 0 0 0'],fmt='%s');

	fp.write('plotcrit ./mckg/cc%d_%d \n \n'%(glno,jj+1));
    fp.close();
     
    print "Running mckg/gllens%d.in on %d"%(glno,num);
    call('%s mckg/gllens%d.in > mckg/gltmpout_%d'%(lenscode,glno,glno),shell=1);
   
    for jj in range(nsrc):
	u1,v1=np.loadtxt('mckg/cc%d_%d'%(glno,jj+1),usecols=(2,3),unpack=True);
	umin=int(min(u1));
	umax=int(max(u1)+1);
	vmin=int(min(v1));
	vmax=int(max(v1)+1);

	call("head -6 mckg/gllens%d.in > mckg/g2lens%d_%d.in"%(glno,glno,jj+1),shell=1);  
	call("head -%d mckg/gllens%d.in | tail -%d >> mckg/g2lens%d_%d.in"%(((jj+1)*nline)+6,glno,nline,glno,jj+1),shell=1);  
        fp=open('mckg/g2lens%d_%d.in' % (glno,jj+1),'a');
        for mm in range(100):	
            sx=np.random.uniform(umin,umax);
            sy=np.random.uniform(vmin,vmax);
            fp.write('findimg %f %f\n'%(sx,sy));
        np.savetxt(fp,['quit'],fmt='%s');
	fp.close();
	 
        print "Running mckg/g2lens%d.in on %d"%(glno,num);
        call('%s mckg/g2lens%d_%d.in > mckg/g2tmpout%d_%d'%(lenscode,glno,jj+1,glno,jj+1),shell=1);
    
    ii=0;
    srcx=np.arange(nsrc)*0.;
    srcy=np.arange(nsrc)*0.;
    smag=np.arange(nsrc)*0.;
    summ=np.arange(nsrc)*0.;
    imno=np.arange(nsrc)*0;
     
    for ii in range(nsrc):
        ip=ii+1;
	imgnos=np.random.uniform(0,1);
	## Loop to select a source that has 2 or more lensed images such that a
	## double is chosen 50% of the times
	cond=1;
	while(cond<3):
	    if(imgnos<0.5):
		call("awk -f scr2img mckg/g2tmpout'%d'_'%d' > mckg/tp_%d"%(glno,ip,glno),shell=True);
		## if no doubly images sources (unlikely to happen), try sources with no. of lensed
		## images >2
		if( os.path.getsize("mckg/tp_%d"%(glno))<=0):
		    imgnos=0.6;
		    cond=cond+1;
		else:
		    rwno=np.loadtxt("mckg/tp_%d"%(glno),unpack=True);
		    cond=10;
	    else:    
		call("awk \'{if($3~/images:/ && $2>2) print NR}\' mckg/g2tmpout'%d'_'%d' > mckg/tp_%d"%(glno,ip,glno),shell=True);
		## if no sources with no. of lensed images >2, try to find
		## doubly imaged sources
		if( os.path.getsize("mckg/tp_%d"%(glno))<=0):
		    imgnos=0.4;
		    cond=cond+1;
		else:
		    rwno=np.loadtxt("mckg/tp_%d"%(glno),unpack=True);
		    cond=10;
        if(cond==3):
	    print "No lensed images found for lens ",ii,"with id",gid,"..continuing..";
	    continue;
        else: 
	    ## More than one eligible sources found 
	    if(rwno.ndim>0):    
		indx=np.random.randint(0,rwno.size);
		rwn=int(rwno[indx]);
		
		filename=('mckg/g2tmpout%d_%d'%(glno,ip));
		print "gllens%d.in"%(glno),filename,rwn;
		
		## Extract source position and no. of lensed images
		call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d \\n\",$2); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s > mckg/tmp_%d"%(filename,glno),shell=1);
		## Extract source position, no. of lensed images, row no. of the
		## first lensed image and lens galaxy id
		call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d %d '+"%d"%(gid)+' \\n\",$2,NR+1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s >> mckg/srcinfo_%d"%(filename,glno),shell=1);
		xx,yy,jj=np.loadtxt('mckg/tmp_%d'%(glno),unpack=True);
		call("awk \'{if(NR>'%d' && NR<='%d') print $3}\' '%s'  >  mckg/tmp2_%d" % (rwn,rwn+jj,filename,glno),shell=True);
		magn=np.loadtxt('mckg/tmp2_%d'%(glno),unpack=True);
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
		
		filename=('mckg/g2tmpout%d_%d'%(glno,ip));
		print "gllens%d.in"%(glno),filename,rwn;
		call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d \\n\",$2); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s > mckg/tmp_%d"%(filename,glno),shell=1);
		call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d %d '+"%d"%(gid)+' \\n\",$2,NR+1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s >> mckg/srcinfo_%d"%(filename,glno),shell=1);
		xx,yy,jj=np.loadtxt('mckg/tmp_%d'%(glno),unpack=True);
		srcx[ii]=xx;
		srcy[ii]=yy;
		imno[ii]=int(jj);

		call("awk \'{if(NR>'%d' && NR<='%d') print $3}\' mckg/g2tmpout'%d'_'%d' >  mckg/tmp2_%d" % (rwn,rwn+jj,glno,ip,glno),shell=True);
		magn=np.loadtxt('mckg/tmp2_%d'%(glno),unpack=True);
		magnsrt=sorted(np.abs(magn),reverse=1);
		magsum=sum(np.abs(magn));
		smag[ii]=magnsrt[1]; 
		summ[ii]=magsum; 
	    
    return srcx,srcy,smag,summ,imno;

        

## Gravlens file to generate simulated galaxy lensed images
def genimg_gg(kappa,rscale,reinst,gposx,gposy,ell,ell_pa,su,sg,sr,si,sz,srcx,srcy,ell_s,ell_pa_s,reff_s,zpt,gid,indxno,nsrc,num):

    ## For various bands
    band=['u','g','r','i','z'];
    
    ## Create the gravlens file for each lens
    fp=open('gout/gglens_%d.in' % (gid),'w');

    np.savetxt(fp,['set rscale=15'],fmt='%s');
    np.savetxt(fp,['set ngrid1=40'],fmt='%s');
    np.savetxt(fp,['set ngrid2=40'],fmt='%s');
    np.savetxt(fp,['set verbose=0'],fmt='%s');
    np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');

    lcno=gposx.size+1;
    fp.write('startup %d 1 \n'%(lcno));
    fp.write('nfw %f 0 0 0 0 0 0 %f 0 0 \n' %(kappa,rscale));
    for mm in range(gposx.size):
	fp.write('alpha %f %f %f %f %f 0 0 0 0 1 \n' %(reinst[mm],gposx[mm],gposy[mm],ell[mm],ell_pa[mm]));
    for mm in range(lcno): 
	np.savetxt(fp,['0 0 0 0 0 0 0 0 0 0'],fmt='%s');

    ## Loop over each band
    for kk in range(5):
	if(kk==0):
	    mag=su;
	elif(kk==1):
	    mag=sg;
	elif(kk==2):
	    mag=sr;
	elif(kk==3):
	    mag=si;
	elif(kk==4):
	    mag=sz;

	flux=10.**(-0.4*(mag-zpt));
	
	## Use deVaucouler's profile, since the size-lum relation is used
	## for this profile and hence, the reff corresponds to the half-light
	## radius for that profile
	np.savetxt(fp,['setsource 1 1'],fmt='%s');
	fp.write('sersic %f %f %f %f %f %f 0 0.5 macro \n' %(flux,srcx,srcy,ell_s,ell_pa_s,reff_s));
	np.savetxt(fp,[' 0 0 0 0 0 0 0 0 '],fmt='%s \n');

	fp.write('SBmap2 -100 100 201 -100 100 201 1 imout_%s_%d_%s.fits 3 \n' % (indxno,gid,band[kk]));
	fp.write('convolve1 imout_%s_%d_%s.fits 3  psfcfh_%s.fits 3 imoutp_%s_%d_%s.fits 3 \n \n' % (indxno,gid,band[kk],band[kk],indxno,gid,band[kk]));

    np.savetxt(fp,['quit'],fmt='%s');
    fp.close();


## Calculate size of the source assuming size-luminosity relation
## using the Bernardi et al. 2003 eqns given in Oguri 2006
def srcsize(mapp,zsrc,pixsc):
    ## mapp is in g-band
    ## Since H0=70 no h70 factors required
    ## distances are converted from h^-1Mpc to Mpc
    Dlum=cc.Dlofz(zsrc)/p.hval;
    Mabs =mapp-5*log10(Dlum)-25;
    Lum_src=10**(-0.4*(Mabs-5.48));
    Da=cc.Daofz(zsrc)*1.e3/p.hval;

    Lrat= Lum_src/10**10.2;
    Reff= 10**0.52*Lrat**(2./3.) * 1./(1+zsrc)**2;
    ## Half-light radius converting from kpc->radian->arcsec->pix
    return (Reff/Da)*(180.*3600/pi/pixsc);
