#!/usr/bin/python
import pyfits
from math import *
import numpy as np
import sys

##Reads finp FITS file with mean counts/sec and uses the exptime from the header. Creates FITS file fout by taking into account the Poisson noise from the detector. 
##Example: ./poisson.py finp fout";
#def poinoise(fname1):
if(len(sys.argv)!=2):
    print "./poisson.py finp.fits";
    sys.exit(0);
else:
    fname1=sys.argv[1];

hdulist=pyfits.open(fname1,mode='update');
##hdulist=pyfits.open(fname1);
hdulist.info();

data=hdulist[0].data;
try:
    exptime=hdulist[0].header['EXPTIME'];
except KeyError:
    print "Running poisson.py on %s"%(fname1);
    exit(0);

data=data*exptime;
data=np.abs(data);
data=np.random.poisson(data)*1.0;
data=data/exptime;
hdulist[0].data=data;
hdulist.flush();
##hdulist.writeto(fname2);
hdulist.close();
