#!/usr/bin/env python
## Custom settings can be made by changing this file. However, 
## see "IMPORTANT NOTE" in the file ../README
import sys
from math import *

## Setting up Cosmology
##################
## Add the full path to aum_mini and the configfile below
ff=open("/disk2/anupreeta/simcode/aum_mini/configfile");
for line in ff:
    linecont=line.strip("\n")[1:];
    sys.path.append("/disk2/anupreeta/simcode/aum_mini/"+linecont);

import cosmology as c

p=c.cosmo();
p.Om0=0.26;
p.Omk=0.0;
p.w0=-1.0;
p.hval=0.72;
p.Omb=0.0469;
p.th=1.0093*2.7;
p.s8=0.8;
p.nspec=0.95;
p.ximax=log10(8.0);
p.cfac=1.0;

cc=c.cosmology(p);

## Some standard paths and constants
##################
## Input foreground and background catalogs
lenscatalog='../../catalogs/test_clus.txt'
#lenscatalog='../../catalogs/inp_gal_frg_cat_w2.txt'
bkggalcatalog="../../catalogs/inp_gal_bkg_cat.txt"

## Catalog with Redshift, M200, Mvir, Rs, Rhos which can be used to calculate
## kappa, the convergence 
fileforkappa="FileFor_Kappa_s_R_s.dat"

## Path to keeton's lensmodel executable
lenscode="/home/anupreeta/soft/lensmodel"

## Constants
gee=4.2994e-9;   ##Mpc/Msun (Km/s)**2
cspd=299792.458; ##km/s

## Pixel scale
pixsc=0.186;

## Zeropoint for the survey
zpt=30;

## No. of processors to use when running this code
Nproc=22;

## See srcprop.py
##################
## Set redshift and magnitude limits and factor by which to boost the source
## number density to get enhanced lens cross-section
zmin_gal=1.0;
zmax_gal=4.0;
mmin_gal=21.0;
mmax_gal=25.5; 
boost_csect_gal=40;

## See main_clus.py
##################
## Set magnitude limit on 2nd brightest image and total flux
## of all the lensed images
mi_lim_cfht_g=23.0; 
mi_totbrt_g=19.; 
    
## Set the range within which to match magnitudes,colors and redshift of real
## galaxies in order to extract colors 
magdiff=5*0.1;

zdiffg_strict=3.;
zdiffg_relaxed=20.;

gicolorg_strict=1.5;
gicolorg_relaxed=1.5;

## Draw uniformly from random ellipticities and PA for the background source
## using the following limits
bkggal_ell_low=0.1;
bkggal_ell_high=0.6;
bkggal_ellpa_low=0.;
bkggal_ellpa_high=180.;

## See getimpos.py
##################
## Keep only those cluster-scale lenses such that their Reinst > Reinst_min
## in units of arcsec, exclude those lenses where the BGG is brighter than
## magi_upperlim and extract randomly Ncl_lens no. of cluster lenses to be
## used in the final catalog
Reinst_min=2.0;
magi_upperlim=15.7;
Ncl_lens=30;
