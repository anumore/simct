#!/usr/bin/env python
import numpy as np
from astropy.io import fits
import random, time
from glob import glob
from math import *
from astropy.table import Table
from subprocess import call
#from mpi4py import MPI
import os
import sys

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


def add_Poisson_noise(filename,kwargs_source,kwargs_deflector,kwargs_survey,outdir,exptime):
    band = kwargs_survey.get('bands')

    file_table = Table.read(filename, format='csv')
    ra,dec,gid = file_table['galaxy_ra'],file_table['galaxy_dec'],file_table['galaxy_id']



    for ii in range(len(ra)):


        subst="%d"%(gid[ii])
        rrh,ddh=deg2hms(ra[ii],dec[ii])
        fid="%s%s"%(rrh,ddh)
        print("##\n For gal ## %d"%(gid[ii]),ra[ii],dec[ii],fid)
        if(os.path.isfile("gout/imoutp_%s_%s.fits"%(subst,band[0])) and not os.path.isfile("%s/imoutf_%s.fits"%(outdir,fid))):
  
            for kk in range(len(band)):
      
                if(kk==0):
                    img_pois=np.zeros((5,101,101), dtype="float32")

                hdulist=fits.open("gout/imoutp_%s_%s.fits"%(subst,band[kk]))
                print('data shape',hdulist[0].data.shape)
                print("###  Adding imoutp_%s_%s to gal_%s_%s"%(subst,band[kk], fid,band[kk]))
                #  ### Add  Poisson noise
                simdata=hdulist[0].data
                simdata=simdata*exptime[kk]
                simdata=np.abs(simdata)
                #np.random.seed(myseed)
                simdata=np.random.poisson(simdata)*1.0
                simdata=simdata/exptime[kk]
            
                img_pois[kk,:,:]=simdata[0:101,0:101]
                print('Final image size: ',img_pois.shape)

            hdu = fits.PrimaryHDU(img_pois)
            hdulist = fits.HDUList([hdu])
            hdulist.writeto("%s/imoutpois_%s.fits"%(outdir,fid),overwrite=True)
            
           #else:
           #    print "%s/imoutf_%s exists already"%(outdir,fid);
        else:
            print("gout/imoutp_%s_%s.fits"%(subst,band[0]),"does not exist or %s/imoutf_%s exists already"%(outdir,fid))


