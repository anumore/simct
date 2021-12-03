#!/bin/bash  
###########################################################
## Primary custom settings: input_qg.py 
###########################################################
## All tasks are divided over Nproc no. of processors, change Nproc in .py files
## here based on the no. of processors available your machine

## Create "Crosssect.dat", if it does not exist
###########################################################
 ./mkcrosssect.py 

## Use foreground galaxy catalogs to extract those galaxies which will act as
## lenses and generate corresponding background (bkg) galaxy or quasar properties
## Also, uses Crosssect.dat and Cosmology.dat
## Extract lenses with Reinst between 1.5 and 2 arcsec only
###########################################################
 rm GAL*.txt QSO*.txt LOG*.txt
 ./main.py

## Nproc no. of output files are created.
## Combine them to make a single catalog for bkg galaxies and for bkg quasars
###########################################################
 cat GAL*txt > all_gal.txt
 cat QSO*txt > all_qso.txt

## Create output directories for galaxies and quasars
###########################################################
## mkdir gout qout

## Copy PSF FITS files to the output dir with galaxy lenses
###########################################################
## cp /path_where_psf_files_were_stored/psfcfh_?.fits gout/

## Create directories for intermediate outputs
###########################################################
 mkdir gal_qso mckq mckg 

 mv G*txt Q*txt L*txt gal_qso/
 cp all_gal.txt all_qso.txt gal_qso/

## Ensure that the same lens galaxy does not get assigned both a background
## galaxy and a quasar 
###########################################################
#awk '{print $1}' all_gal.txt | sort -u -n > uniqidg
#awk '{print $1}' all_qso.txt | sort -u -n > uniqidq
#./selqso.py
# 
## Generate an extended lens catalog with all parameters necessary for
##generating simulated lensed images 
###########################################################
#./mkcatalog.py

## Extract a single quasar and all its parameter values for a given lens and run
## Gravlens to generate lensed quasar images and a final catalog with lens+qso
## properties
###########################################################
#rm mckq/* fpar1_*
#rm qout/gq*in qout/imout*fits
#./main_qso.py
#cat fpar1_?.txt fpar1_??.txt | sort -n -k1 -u > finalpar1.txt
# 
## Extract a single galaxy and all its parameter values for a given lens and 
## create Gravlens input file to generate lensed galaxy images and a final
## catalog with lens+galaxy properties
###########################################################
#rm mckg/* fpar0_*
#rm gout/gg*in gout/imout*fits
#./main_gal.py
#cat fpar0_?.txt fpar0_??.txt | sort -n -k1 -u > finalpar0.txt
#mv fpar?_*txt gal_qso/

## Run Gravlens to generate lensed images using .py file below
## Need Keeton's code along with the path to the executable
###########################################################
#cd gout
## Create link if it does not exist
#ln -s ../rungl.py .
#./rungl.py
#cd ..
 
## Note all galaxy (or quasar) lensed images are stored in gout (or qout) including the
## input files on which Gravlens is run

## Exclude those lens galaxies which are labelled as BCGs in the cluster catalog
###########################################################
#./remcllens.py

###Clean all intermediate and final output files to start afresh
###########################################################
##rm *.txt uniq*  qout/* gout/g*in gout/imout*fits gout/LO*dat
##rm -r gal_qso mckq mckg 
