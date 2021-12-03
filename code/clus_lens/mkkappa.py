#!/usr/bin/env python
import numpy as np
from math import *
from input_g import *

try:
    with open('FileFor_Kappa_s_R_s.dat'):
	print "FileFor_Kappa_s_R_s.dat already exists, exiting without creating a new file. If cosmological model has been modified, please remove and recreate this file."
	exit();

except IOError:
    print "File FileFor_Kappa_s_R_s.dat does not exist, creating it now";
    fp=open("FileFor_Kappa_s_R_s.dat","w");
    for ii in range(1500):
	zz=0.0+1.5*(ii)/(1499.);
	for jj in range(100):
	    xM=9.0+jj*7./99.;
	    m200=10**xM;
	    Mvir=c.dp();
	    Rvir=c.dp();
	    cvir=c.dp();
	    cc.modelNFWhalo_com(m200,zz,Mvir,Rvir,cvir);
	    cv=cvir.value();
	    rs=Rvir.value()/cv;
	    rhos=Mvir.value()/(4*pi*rs**3.)/(log(1+cv)-cv/(1+cv));
	    fp.write('%3.2e %7.6e %7.6e %6.5e %6.5e \n'%(zz,log10(m200),log10(Mvir.value()),log10(rs),log10(rhos)));

    fp.close();
