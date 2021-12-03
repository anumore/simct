#!/usr/bin/env python
import numpy as np
from math import *
from input_qg import *
import srclensprop as slp

try:
    with open('Crosssect.dat'):
	print "Crosssect.dat already exists, exiting without creating a new file"
	exit();
except IOError:
      print "File Crosssect.dat does not exist, creating it now";
      qi=np.arange(0.1,1.0,0.001);
      csecti=qi*0.0;
      for ii in range(qi.size):
	csecti[ii]=slp.getcrosssect(slp.fid_b_I/sqrt(qi[ii]),qi[ii]);
      np.savetxt("Crosssect.dat",np.transpose([qi,csecti]));

