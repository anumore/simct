#!/usr/bin/env python
## Custom settings can be made by changing this file. However, 
## see IMPORTANT NOTE in the file ../README
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
lenscatalog='../../catalogs/test_gal.txt'
##lenscatalog='../../catalogs/inp_gal_frg_cat_w2.txt'
bkgqsocatalog="../../catalogs/inp_qso_bkg_cat.txt"
bkggalcatalog="../../catalogs/inp_gal_bkg_cat.txt"

## Path to keeton's lensmodel executable
lenscode="/home/anupreeta/soft/lensmodel"

## Constants
gee=4.2994e-9;
cspd=299792.458;

## Pixel scale
pixsc=0.186;

## Zeropoint for the survey
zpt=30;

## No. of processors to run this code on
Nproc=22;

## See main.py
##################
## Set the min max range for reinst(arcsec) of quasar and galaxy lenses
reinst_ql=0.5
reinst_qh=3.0
reinst_gl=0.7
reinst_gh=4.0
zlens_min=0.1
zlens_max=1.1

## See srcprop.py
##################
## Set redshift and magnitude limits and factor by which to boost the source
## number density to get enhanced lens cross-section
zmin_qso=1.0;
zmax_qso=5.9;
mmin_qso=21.0;
mmax_qso=25.5; 
boost_csect_qso=1200;

zmin_gal=1.0;
zmax_gal=4.0;
mmin_gal=21.0;
mmax_gal=25.5; 
boost_csect_gal=100;

## See mkcatalog.py
##################
## Limits on the shear strength and PA. Lens models have shear drawn from this
## range
shear_strength_low=0.001;
shear_strength_high=0.02;
shear_pa_low=0.;
shear_pa_high=180.;

## See main_gal.py and main_qso.py
##################
## Set magnitude limit on 2nd brightest image and total flux
## of all the lensed images
mi_lim_cfht_g=23.0; 
mi_totbrt_g=19.; 
mi_lim_cfht_q=23.0; 
mi_totbrt_q=20.; 
    
## Set the range within which to match magnitudes,colors and redshift of real
## background quasars /galaxies in order to extract colors 
magdiff=5*0.1;
zdiffq=0.1;
zdiffg=0.5;

zdiffq_strict=3.*zdiffq;
zdiffq_relaxed=40.*zdiffq;
zdiffg_strict=3.;
zdiffg_relaxed=20.;

gicolorg_strict=1.5;
gicolorg_relaxed=1.5;
gicolorq_relaxed=1.5;

## Draw uniformly from random ellipticities and PA for the background source
## using the following limits
bkggal_ell_low=0.1;
bkggal_ell_high=0.6;
bkggal_ellpa_low=0.;
bkggal_ellpa_high=180.;

## See genimg.py
##################
## FWHM of seeing in arcsec
sig=[0.85,0.78,0.71,0.64,0.68];

## Set image size of the sims
imsize=101;
band=['g','r','i']
exptime=np.array([600.,600.,600.])##,600.])

## See remcllens.py
##################
## Path to the final cluster lens catalog, to exclude lensing galaxies that are 
## common to the quasar lenses or galaxy-scale lenses and cluster lens catalog 
cluslenscatalog="../clus_lens/finpar_g2.txt"
N_lens=30
