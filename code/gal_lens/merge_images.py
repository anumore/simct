#!/usr/bin/env python
import numpy as np
from astropy.io import fits
import random, time
from glob import glob
from math import *
from subprocess import call
#from mpi4py import MPI
import os
import sys
from astropy.table import Table

def deg2hms(rr,dd):
    flag=0
    rr1=floor(rr/15.)
    rr2=floor((rr/15.-rr1)*60)
    rr3=(((rr/15.-rr1)*60) - rr2)*60

    if(dd<0):
        flag=1
        dd=-1*dd

    dd1=floor(dd)
    dd2=floor((dd-dd1)*60)
    dd3=(((dd-dd1)*60)-dd2)*60

    if(flag):
        return "%02d%02d%02d"%(rr1,rr2,floor(rr3)), "-%02d%02d%02d"%(dd1,dd2,floor(dd3))
    else:
        return "%02d%02d%02d"%(rr1,rr2,floor(rr3)), "+%02d%02d%02d"%(dd1,dd2,floor(dd3))

#band1=['G','R','I','Z']
#band=['g','r','i','z']
#myseed = 109024

def addto_realimage(filename,kwargs_source,kwargs_deflector,kwargs_survey,fitsdir,outdir,exptime):
    band = kwargs_survey.get('bands')

    file_table = Table.read(filename, format='csv')
    ra,dec = file_table['galaxy_ra'], file_table['galaxy_dec']
    gid = file_table['galaxy_id']

    for ii in range(len(ra)):

        subst="%d"%(gid[ii])
        rrh,ddh=deg2hms(ra[ii],dec[ii])
        fid="%s%s"%(rrh,ddh)
        print("##\n For gal ## %d"%(gid[ii]),ra[ii],dec[ii])

        hdulist=fits.open("%s/imoutpois_%s.fits"%(outdir,fid))
        print('poisson image shape: ',hdulist[0].data.shape)
        fin_im = hdulist[0].data
        #raw_im = np.zeros_like(fin_im)
        if(os.path.isfile("%s/imoutpois_%s.fits"%(outdir,fid)) and not os.path.isfile("%s/imoutf_%s.fits"%(outdir,fid))):

            for kk in range(len(band)):
                try:
                    hdulist1=fits.open("%s/hsc-%s/J%s_%s.fits"%(fitsdir,band[kk],fid,band[kk]))
                    #hdulist1=fits.open("%s/hsc-%s/gal_%s_%s.fits"%(fitsdir,band[kk],fid,band[kk]))

                    hsc_galim=hdulist1[1].data

                    axsize1=hdulist1[1].header['NAXIS1']
                    axsize2=hdulist1[1].header['NAXIS2']       
                    print(axsize1,axsize2)

                    if(axsize1==120 and axsize2==120):
                        print("###  Adding imoutp_%s_%s to gal_%s_%s"%(subst,band[kk], fid,band[kk]))             
                        fin_im[kk,:,:]=fin_im[kk,:,:]+hsc_galim[10:-9,10:-9]
                        #raw_im[kk,:,:] = hsc_galim[10:-9,10:-9]
                        print('Final image size: ',fin_im.shape)
                        flag = 1
                    else:
                        print("## gal_%s is small, skipping "%(fid))
                        flag = 0
                        continue

                except IOError:
                    print("%s/hsc-%s/J%s_%s.fits"%(fitsdir,band[kk],fid,band[kk]),"does not exist")
                    continue
            if(flag==1):
                hdu = fits.PrimaryHDU(fin_im)
                hdulist = fits.HDUList([hdu])
                hdulist.writeto("%s/imoutf_%s.fits"%(outdir,fid))
                #fits.writeto("%s/im_raw_%s.fits"%(outdir,fid),raw_im)
            else:
                continue
           #else:
           #    print "%s/imoutf_%s exists already"%(outdir,fid);
        else:
            print("%s/imoutpois_%s.fits"%(outdir,fid),"does not exist or %s/imoutf_%s exists already"%(outdir,fid))


