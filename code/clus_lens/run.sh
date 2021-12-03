#!/bin/bash  
###########################################################
## Primary custom settings: input_qg.py 
###########################################################
## All tasks are divided over Nproc no. of processors, change Nproc in .py files
## here based on the no. of processors available your machine

## Create "FileFor_Kappa_s_R_s.dat", if it does not exist
###########################################################
 ./mkkappa.py > /dev/null 2>&1 

## Use cluster catalogs to extract those clusters which will act as
## lenses and generate corresponding background (bkg) galaxy properties
## Also, uses Cosmology.dat
###########################################################
 rm GAL*.txt LOG*.txt mck1/*
 ./main.py

## Nproc no. of output files are created.
## Create a single catalog for bkg galaxies per each lens and create a single
## catalog with info on BCG+members   
###########################################################
 cat GAL*txt > all_gal.txt
 cat LEN*txt > lenscat.txt

## Create directories for intermediate outputs
###########################################################
 mkdir galt mck1 mckg 

 mv G*txt L*txt galt/
 cp all_gal.txt lenscat.txt galt/
  
## Generate an extended lens catalog with all parameters necessary for
## generating simulated lensed images 
###########################################################
#./mkcatalog.py
#  
## Extract a single bkg galaxy and all its parameter values for a given lens and 
## create Gravlens input file to generate lensed galaxy images and a final
## catalog with lens+galaxy properties
###########################################################
#rm mckg/* 
#rm gout/gg*in gout/imout*fits
#./main_clus.py
#cat fpar_?.txt fpar_??.txt | sort -n -k1 -u > finalpar0.txt
#mv fpar_* galt/

## Run Gravlens to generate lensed images using .py file below
## Need Keeton's code along with the path to the executable
###########################################################
#cd gout
#ln -s ../rungl.py .
#./rungl.py
#cd ..

## Get source info for each lens
###########################################################
#awk 'function abs(x){return x<0?-x:x}{printf("%e %e %d \n",abs($14),abs($15),$1)}' finalpar0.txt > sinf0
#rm sinfgrepout
#awk '{ print "grep \""$1" ."$2"\" mckg/g2tmpout*  >> sinfgrepout"}' sinf0 > rgreplist
#sh rgreplist
#awk -F ":" '{if($2~/# source/) print $1,$2}' sinfgrepout > sg2list
#paste sg2list sinf0 | awk '{print $1,$2,$3,$6,$7,$8}' > srcinf_inp
#cat mckg/srcinfo_* > srcinfoall
#rm sinf0 rgreplist sinfgrepout sg2list

## Extract image positions to calculate Reinst (Einstein radius)
## Keep those lenses for which Reinst>2.0 arcsec
###########################################################
#./getimpos.py


###Clean all intermediate and final output files to start afresh
###########################################################
##rm *.txt srcinf* img* gout/gg*in gout/imout*fits
##rm -r galt mck1 mckg  

