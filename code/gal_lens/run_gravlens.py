#!/usr/bin/env python
from math import *
import numpy as np
import sys
#from StringIO import StringIO
from subprocess import call
#import multiprocessing
#from mpi4py import MPI
import os
import time
import glob
from astropy.table import Table
#sys.path.append("../")

###########################################################
## PURPOSE:
## Run lensmodel on gout/*.in with sources as bkg galaxies 
###########################################################

def rungl(lenscode,filename):

    """
    Function to execute lensing simulations for galaxies.

    Parameters:
    lenscode (str): path to the 'lesncode' software to run the lensing simulation.
    filename (str): Input filename to run lensmodel.
    """

    param_table = Table.read(filename, format='csv')
    galaxy_id = param_table['galaxy_id']


    print(len(galaxy_id),' galaxies', len(np.unique(galaxy_id)),' unique')
    for ii in range(len(galaxy_id)):

        subst = f"{galaxy_id[ii]}"
        
        # Check if the output file already exists
        output_file = f"gout/imoutp_{subst}_g.fits"
        input_file = f"gout/gglens{subst}_im.in"
        log_file = f"LOG.{subst}.dat"


        if(not os.path.isfile(output_file)):
            if(os.path.isfile(input_file)):
                # If input file exists, run the lensing simulation
                print(f"{lenscode} {input_file} > {log_file}")
                call(f"{lenscode} {input_file} > {log_file} 2>&1", shell=True)
            else:
                # If input file doesn't exist, log a message
                print(f"{input_file} doesn't exist")
        else:
            # If output file already exists, log a message
            print(f"{output_file} exists")


