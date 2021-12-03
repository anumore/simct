#!/usr/bin/env python
from math import *
import numpy as np
import sys
from subprocess import call

fldid,indxno=np.loadtxt('fieldid_sort',dtype={'names': ('fldid','indxno'), 'formats': ('S13', 'S3',)}, usecols=(1,0),unpack=True);

for ii in range(fldid.size):         
    for band in ['u','g','r','i','z']:
	command="cphead fitsfiles/CFHTLS_W_%s_%s_T0007_MEDIAN.fits gout/imoutp_%s_*_%s.fits exptime"%(band,fldid[ii],indxno[ii],band);
	print command;
	call(command,shell=1);
