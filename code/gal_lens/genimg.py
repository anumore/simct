#!/usr/bin/env python
from math import *
import numpy as np
from subprocess import call
import os
from astropy.io import fits
#from cosmo_params import *
from scipy import interpolate

## Spatial scales are in pixels and flux is converted to counts/sec
## imsize is the size of each simulated image
## grizY are the 5 bands that are being simulated here

## Gravlens file to generate simulated quasar lensed images

#imsize = 101
#lenscode = "/global/homes/v/vibhore/soft/lensmodel"
#psfpath       = "../inppsf"

def genimg_gq(reinst,ell,ell_pa,sh_str,sh_pa,sg,sr,si,sz,sy,srcx,srcy,zpt,gid,num,lenscode,psfpath,imsize):

    band=['g','r','i','z','y']
    ## for various bands
    mag=[sg,sr,si,sz,sy]

    ## create the gravlens file for each lens
    fvar=''
    fvar += ('qout/gqlens%s_im.in' % (gid))
    #fvar=('gout/gl.in');
    fp=open(fvar,'w')

   #np.savetxt(fp,['set rscale=15'],fmt='%s');
   #np.savetxt(fp,['set ngrid1=40'],fmt='%s');
   #np.savetxt(fp,['set ngrid2=40'],fmt='%s');
   #np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');
    fp.write("set rscale=15 \n set ngrid1=40 \n set ngrid2=40 \n set omitcore=0.02 \n")

   #np.savetxt(fp,['startup 1 1'],fmt='%s');
    fp.write("startup 1 1\n")
    fp.write('alpha %f 0 0 %f %f %f %f 0 0 1 \n' %(reinst,ell,ell_pa,sh_str,sh_pa))
   #np.savetxt(fp,[' 0 0 0 0 0 0 0 0 0 0 '],fmt='%s \n');
    fp.write(" 0 0 0 0 0 0 0 0 0 0 \n")
    fp.write('findimg %f %f \n'%(srcx,srcy))
    fp.write("quit")
    fp.close()

    ## Run lensmodel to generate the lensed images
    cmd_str=("%s qout/gqlens%s_im.in | tail -20 > qout/chkout_%d"%(lenscode,gid,num));
    os.system(cmd_str)

    ## Extract the positions and magnifications of the lensed images
    call("awk \'{if($NF~/tdel/) print $2}\' qout/chkout_%d > qout/var_%d"%(num,num),shell=True);
    fo=open('qout/var_%d'%(num),'r')
    imno=int(fo.read()[:-1])
    fo.close()
    call("tail -%d qout/chkout_%d | head -%d > qout/chkout1_%d"%(imno+2,num,imno,num),shell=True);
    imx,imy,immag=np.loadtxt('qout/chkout1_%d'%(num),usecols=(0,1,2),unpack=True);

    flx=np.arange(len(band))*0
    psfdict={}
    im_arr={}

    imx=(imx)+(imsize/2)
    imy=(imy)+(imsize/2)
    for kk in range(len(band)):
        ## Convert the magnitude to flux in counts
        flx[kk]=10.**(-0.4*(mag[kk]-zpt))
        
        ## Read the psf and interpolate over all pixels
       #hdulist=pyfits.open("DES2145+0001_3019196822_psf_g.fits")
       #hdulist=pyfits.open("%s/%s/%s_%s_psf_%s.fits"%(psfpath,gfld,gfld,gid,band[kk]))
        hdulist=fits.open("%s/psf_hsc_%s.fits"%(psfpath,band[kk]))
       #hdulist=pyfits.open("../inppsf/psf_i.fits")
        scidata = hdulist[0].data
        psfimsize=hdulist[0].header['naxis1']
        ay=np.arange(psfimsize)-(psfimsize/2)
        ax=np.arange(psfimsize)-(psfimsize/2)
        func=interpolate.interp2d(ax, ay, scidata.flatten(), kind='cubic')
        psfdict['func_%s'%(band[kk])]=func
       
        ## Initialize array to generate the output FITS image
        arr=np.zeros((imsize,imsize), dtype="float32")
        im_arr['im_arr_%s'%(band[kk])]=arr


   
    ## Create the simulated lensed qso image
    ## For each pixel in the output image, assign the flux that is expected to arise due to magnification of each of the lensed image and psf convolution

    value=np.arange(len(band))
    for nn in range(imx.size):
        xx=imx[nn] 
        yy=imy[nn]
        for ii in range(imsize):
            for jj in range(imsize):
                for kk in range(len(band)):
                    try:
                        value[kk] = flx[kk]*abs(immag[nn])*psfdict["func_%s"%(band[kk])](jj-yy,ii-xx)
                    except:
                        value[kk] = 0.0
                    im_arr['im_arr_%s'%(band[kk])][jj,ii]=im_arr['im_arr_%s'%(band[kk])][jj,ii] + value[kk]

    
    for kk in range(len(band)):
        hdu = fits.PrimaryHDU(im_arr['im_arr_%s'%(band[kk])])
        hdulist = fits.HDUList([hdu])
        hdulist.writeto('qout/imoutp_%s_%s.fits'%(gid,band[kk]))

    ## Write a catalog with the lens ID, Reinst, image multiplicity, image positions and image magnitudes for each band
    fp1=open('qsocatf_%d.txt'%(num),'w')
    fp1.write("%s %f %d"%(gid,reinst*0.186,imno))
    for nn in range(imx.size):
        fp1.write(" %f %f"%(xx,yy))
        for kk in range(len(band)):
            fp1.write(" %f"%(-2.5*log10(np.abs(immag[nn]))+mag[kk]))
            if(nn==imx.size-1 and kk==len(band)-1):
               fp1.write(" \n")
    fp1.close()



## Gravlens file to generate simulated galaxy lensed images
def genimg_gg(reinst,ell,ell_pa,sh_str,sh_pa,sg,sr,si,sz,sy,srcx,srcy,ell_s,ell_pa_s,reff_s,zpt,gid,lenscode,psfpath,imsize):

    ## For various bands
    band=['g','r','i','z', 'y']
    mag=[sg,sr,si,sz, sy]
    ## Create the gravlens file for each lens
    fvar=''
    fvar += ('gout/gglens%s_im.in' % (gid))
    fp=open(fvar,'w')

   #np.savetxt(fp,['set rscale=15'],fmt='%s');
   #np.savetxt(fp,['set ngrid1=40'],fmt='%s');
   #np.savetxt(fp,['set ngrid2=40'],fmt='%s');
   #np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');
    fp.write("set rscale=15 \n set ngrid1=40 \n set ngrid2=40 \n set omitcore=0.02 \n")

    fp.write("startup 1 1\n")
    fp.write('alpha %f 0 0 %f %f %f %f 0 0 1 \n' %(reinst,ell,ell_pa,sh_str,sh_pa))
   #np.savetxt(fp,[' 0 0 0 0 0 0 0 0 0 0 '],fmt='%s \n');
    fp.write(" 0 0 0 0 0 0 0 0 0 0 \n")

    flux=np.arange(len(band))*0
    flux=10.**(-0.4*(np.asarray(mag)-zpt))

    ## Loop over each band
    for kk in range(len(band)):
        
        ## Using deVaucouler's profile, since the size-lum relation is used
        ## for this profile and hence, the reff corresponds to the half-light
        ## radius for that profile
        ## XXX check the size, profile values
        fp.write("setsource 1 1 \n")
        fp.write('sersic %f %f %f %f %f %f 0 1.0 macro \n' %(flux[kk],srcx,srcy,ell_s,ell_pa_s,reff_s))
        fp.write(" 0 0 0 0 0 0 0 0 \n")
        fp.write('SBmap2 -%d %d %d -%d %d %d 1 gout/imout_%s_%s.fits 3 \n' % (imsize/2,imsize/2,imsize,imsize/2,imsize/2,imsize,gid,band[kk]))
        fp.write('convolve1 gout/imout_%s_%s.fits 3  %s/psf_hsc_%s.fits 3 gout/imoutp_%s_%s.fits 3 \n \n' % (gid,band[kk],psfpath,band[kk],gid,band[kk]))
       #fp.write('convolve1 gout/imout_%s_%s_%s.fits 3  %s/%s/%s_%s_psf_%s.fits 3 gout/imoutp_%s_%s_%s.fits 3 \n \n' % (gfld,gid,band[kk],psfpath,gfld,gfld,gid,band[kk],gfld,gid,band[kk]));

   # np.savetxt(fp,['quit'],fmt='%s');
    fp.write("quit")
    fp.close()


