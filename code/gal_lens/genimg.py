#!/usr/bin/env python
from math import *
import numpy as np
from subprocess import call
import os
import pyfits
from input_qg import *

## Spatial scales are in pixels and flux is converted to counts/sec
## 201 pix is the image size of each simulated image
## ugriz are the 5 bands that are being simulated here

## Gravlens file to generate simulated quasar lensed images
def genimg_gq(reinst,ell,ell_pa,sh_str,sh_pa,su,sg,sr,si,sz,srcx,srcy,zpt,gid,indxno,num):

    ## For various bands
    band=['u','g','r','i','z'];

    ## FWHM of seeing in arcsec
  # sig=[0.85,0.78,0.71,0.64,0.68];
    for kk in range(len(band)):
        sig[kk]=sig[kk]/pixsc/(2*sqrt(2*log(2)));
    
    ## Create the gravlens file for each lens
    fvar='';
    fvar += ('qout/gqlens%s_im.in' % (gid));
    fp=open(fvar,'w');

    np.savetxt(fp,['set rscale=15'],fmt='%s');
    np.savetxt(fp,['set ngrid1=40'],fmt='%s');
    np.savetxt(fp,['set ngrid2=40'],fmt='%s');
    np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');

    np.savetxt(fp,['startup 1 1'],fmt='%s');
    fp.write('alpha %f 0 0 %f %f %f %f 0 0 1 \n' %(reinst,ell,ell_pa,sh_str,sh_pa));
    np.savetxt(fp,[' 0 0 0 0 0 0 0 0 0 0 '],fmt='%s \n');
    fp.write('findimg %f %f \n'%(srcx,srcy));	
    np.savetxt(fp,['quit'],fmt='%s');
    fp.close();

    cmd_str=("%s qout/gqlens%s_im.in | tail -20 > qout/chkout_%d"%(lenscode,gid,num));
    os.system(cmd_str);
    call("awk \'{if($NF~/tdel/) print $2}\' qout/chkout_%d > qout/var_%d"%(num,num),shell=True);
    fo=open('qout/var_%d'%(num),'r');
    avar=int(fo.read());
    fo.close();
    ## Extract the lensed image positions and magnifications
    call("tail -%d qout/chkout_%d | head -%d > qout/chkout1_%d"%(avar+2,num,avar,num),shell=True);
    imx,imy,immag=np.loadtxt('qout/chkout1_%d'%(num),usecols=(0,1,2),unpack=True);
    
    imx=imx+100;
    imy=imy+100;

    flxu=10.**(-0.4*(su-zpt));
    flxg=10.**(-0.4*(sg-zpt));
    flxr=10.**(-0.4*(sr-zpt));
    flxi=10.**(-0.4*(si-zpt));
    flxz=10.**(-0.4*(sz-zpt));

    arru=np.zeros((201,201), dtype="float32");
    arrg=np.zeros((201,201), dtype="float32");
    arrr=np.zeros((201,201), dtype="float32");
    arri=np.zeros((201,201), dtype="float32");
    arrz=np.zeros((201,201), dtype="float32");

    ## Calculate flux of the lensed images in each band
    for aa in range(201):
        for bb in range(201):
	    for nn in range(imx.size):
		arru[aa,bb]=arru[aa,bb]+flxu*abs(immag[nn])/(2.*pi*sig[0]**2.) * exp( -((aa-imy[nn])**2.+(bb-imx[nn])**2.)/(2*sig[0]**2.) ) ;
		arrg[aa,bb]=arrg[aa,bb]+flxg*abs(immag[nn])/(2.*pi*sig[1]**2.) * exp( -((aa-imy[nn])**2.+(bb-imx[nn])**2.)/(2*sig[1]**2.) ) ;
		arrr[aa,bb]=arrr[aa,bb]+flxr*abs(immag[nn])/(2.*pi*sig[2]**2.) * exp( -((aa-imy[nn])**2.+(bb-imx[nn])**2.)/(2*sig[2]**2.) ) ;
		arri[aa,bb]=arri[aa,bb]+flxi*abs(immag[nn])/(2.*pi*sig[3]**2.) * exp( -((aa-imy[nn])**2.+(bb-imx[nn])**2.)/(2*sig[3]**2.) ) ;
		arrz[aa,bb]=arrz[aa,bb]+flxz*abs(immag[nn])/(2.*pi*sig[4]**2.) * exp( -((aa-imy[nn])**2.+(bb-imx[nn])**2.)/(2*sig[4]**2.) ) ;

    ## Write FITS images with quasar lensed images as Gaussians of the size of
    ## the seeing 
    print "qout/gqlens%s_im.in qout/imoutp_%s_%d_u.fits \n" % (gid,indxno,gid);
    hdu = pyfits.PrimaryHDU(arru);
    hdulist = pyfits.HDUList([hdu]);
    hdulist.writeto('qout/imoutp_%s_%d_u.fits'%(indxno,gid));

    hdu = pyfits.PrimaryHDU(arrg);
    hdulist = pyfits.HDUList([hdu]);
    hdulist.writeto('qout/imoutp_%s_%d_g.fits'%(indxno,gid));

    hdu = pyfits.PrimaryHDU(arrr);
    hdulist = pyfits.HDUList([hdu]);
    hdulist.writeto('qout/imoutp_%s_%d_r.fits'%(indxno,gid));
    
    hdu = pyfits.PrimaryHDU(arri);
    hdulist = pyfits.HDUList([hdu]);
    hdulist.writeto('qout/imoutp_%s_%d_i.fits'%(indxno,gid));
    
    hdu = pyfits.PrimaryHDU(arrz);
    hdulist = pyfits.HDUList([hdu]);
    hdulist.writeto('qout/imoutp_%s_%d_z.fits'%(indxno,gid));


## Gravlens file to generate simulated galaxy lensed images
def genimg_gg(reinst,ell,ell_pa,sh_str,sh_pa,su,sg,sr,si,sz,srcx,srcy,ell_s,ell_pa_s,reff_s,zpt,gid,indxno,num):

    ## For various bands
    band=['u','g','r','i','z'];
    ## Create the gravlens file for each lens
    fvar='';
    fvar += ('gout/gglens%s_%s_im.in' % (gid,indxno));
    fp=open(fvar,'w');

    np.savetxt(fp,['set rscale=15'],fmt='%s');
    np.savetxt(fp,['set ngrid1=40'],fmt='%s');
    np.savetxt(fp,['set ngrid2=40'],fmt='%s');
    np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');

    np.savetxt(fp,['startup 1 1'],fmt='%s');
    fp.write('alpha %f 0 0 %f %f %f %f 0 0 1 \n' %(reinst,ell,ell_pa,sh_str,sh_pa));
    np.savetxt(fp,[' 0 0 0 0 0 0 0 0 0 0 '],fmt='%s \n');

    ## Loop over each band
    for kk in range(len(band)):
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
	
	## Using deVaucouler's profile, since the size-lum relation is used
	## for this profile and hence, the reff corresponds to the half-light
	## radius for that profile
	np.savetxt(fp,['setsource 1 1'],fmt='%s');
	fp.write('sersic %f %f %f %f %f %f 0 0.5 macro \n' %(flux,srcx,srcy,ell_s,ell_pa_s,reff_s));
	np.savetxt(fp,[' 0 0 0 0 0 0 0 0 '],fmt='%s \n');

	fp.write('SBmap2 -100 100 201 -100 100 201 1 imout_%s_%d_%s.fits 3 \n' % (indxno,gid,band[kk]));
	fp.write('convolve1 imout_%s_%d_%s.fits 3  psfcfh_%s.fits 3 imoutp_%s_%d_%s.fits 3 \n \n' % (indxno,gid,band[kk],band[kk],indxno,gid,band[kk]));


    np.savetxt(fp,['quit'],fmt='%s');
    fp.close();

