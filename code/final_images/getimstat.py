#!/usr/bin/python2.7
import pyfits
from math import *
import numpy as np
import sys

### Add min, max, median and mean of the image to its header
if(len(sys.argv)!=2):
    print "./getimstat.py finp.fits";
    sys.exit(0);
else:
    fname1=sys.argv[1];

###def addstats(filefits):
hdulist=pyfits.open(fname1,mode='update');
data=hdulist[0].data;

#### exclude pixels with value 0, when estimating the following stats
arr=hdulist[0].data;
data=arr[arr!=0];

mean=np.mean(data);
median=np.median(data);
minim=np.min(data);
maxim=np.max(data);

hdulist[0].header['MIN']=minim;
hdulist[0].header['MAX']=maxim;
hdulist[0].header['MEDIAN']=median;
hdulist[0].header['MEAN']=mean;

hdulist.close();
