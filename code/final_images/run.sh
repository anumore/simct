#!/bin/bash
 
function partone() {
    fil1=$1
    dirfitsname=$2
    cp $fil1 ./finalpar.txt
    ##############################

    sort -k5 finalpar.txt > finalpar_srt.txt
    ./genpixlist.py
    rm fpxlst
    sh rgetpix
    paste inpgal0 fpxlst > inpgal.txt
    awk '{printf(" %s %f %f %s %s %f %f %04d \n",$1,$2,$3,$4,$5,$6,$7,50*int(($6-1)/386)+int(($7-1)/386)+1)}' inpgal.txt > pximtag
    ./gencombim2.py $dirfitsname inpgal.txt > logcomb

    ## Run the following two commands either here or right at the beginning
    ## Copy exptime from the parent tiles to their respective sim images 
    # sh rchd1
    ## Add poisson noise to the simulated arc pixels (gout/imoutp*) before merging with the real data in the parent tile
    #./runpoi.py > poiout 
    
    ### Check if files look alright 
    echo "imdir/idpxlst"
    awk '{print $2}' imdir/idpxlst | sort -u | wc -l
    echo "finalpar.txt"
    awk '{print $1}' finalpar.txt | sort -u | wc -l
    echo "inpgal.txt"
    awk '{print $1}' inpgal.txt | sort -u | wc -l
    wc -l logcomb

    cd imdir
     ./runcombim.py
    cd ..
    
    ### Check if files look alright 
    echo "imdir/idpxlst"
    awk '{print $1}' imdir/idpxlst | sort -u | wc -l
    echo "blanksims/*_i.fits ", `ls -l blanksims/*_i.fits | wc -l`
    echo "outfits1/*_i.fits ", `ls -l outfits1/*_i.fits | wc -l`
}

function parttwo(){
    outdirname=$1
    mkdir $outdirname
    ./mkpng_blnk.py
    ./mkpng_cfht.py
    ./mkpng_mask.py
    ./mkpng_crsh.py $outdirname/png
    echo "blanksims/*png ",`ls -l blanksims/*png | wc -l`
    echo "outfits1/*_o_gri.png ",`ls -l outfits1/*_o_gri.png | wc -l`
    echo $outdirname"/png/*png ",`ls -l $outdirname/png/*png | wc -l`
}

function partthree(){
    outdirname=$1
    rm -r $outdirname/blanksims $outdirname/outfits1
    cp finalpar.txt inpgal.txt imdir/idpxlst $outdirname/ 
    cp -r outfits1 $outdirname/
    cp -r blanksims $outdirname/
    ls -l $outdirname/blanksims/*_g.fits | wc -l
    ls -l $outdirname/outfits1/*_g.fits | wc -l
    ls -l $outdirname/outfits1/*_o_gri.png | wc -l
}

#############################################################################
## DIR and FILES NEEDED 
## fitsfiles outfits gout blanksims imdir  outfits1
## finalpar.txt: final lens catalog
## cutoutlist: xmin,xmax,ymin,ymax,cutout_no - pixel values of each contiguous cutout extracted from the parent survey tile
## gout: simulated lens images are stored here
## fitsfiles: the original survey tiles (FITS) are stored here
## outfits: uniform size small cutouts (FITS) of contiguous survey region are stored
## here locally
#############################################################################

###############################
## RUN THIS FOR THE FIRST TIME TO SET UP ALL THE INPUT FILES
###############################

# mkdir fitsfiles outfits gout  blanksims imdir  outfits1

## Create a link to the color composite making code - HumVI/compose.py
# ln -s /path_where_compose.py_from_HumVI_is_located/compose.py .

## Create small cutouts of survey from the parent tiles and creates
## the file "cutoutlist"
## NOTE: this is the same dir path that will be provided to partone()
# ./mkcutouts.py /path_of_remote_dir_where_real_image_cutouts_will_be_stored 440 386 19354

## Create local copies of the simulated images and lens catalogs
## E.g. for galaxy-galaxy lens
#cp ../gal_lens/gout/imoutp*fits gout/
#cp ../gal_lens/finalpar_g.txt ./finalpar.txt
  
## Run the following two commands either here or from within partone() 
## 1. Copy exptime from the parent tiles to their respective sim images 
#./copyhdr.py
## 2. Add poisson noise to the simulated arc pixels (gout/imoutp*) before merging with the real data in the parent tile
#./runpoi.py > poiout 

###############################
## MAIN
###############################
## Input param file and output dir name
#fil1=./finalpar.txt
#dirfitsname=/path_of_remote_dir_where_real_image_cutouts_will_be_stored
#outdirname=galw1
#
### Generates and run scripts which allow you to merge the simulated lensed images with the real
### survey image cutouts (FITS files) using the lensing galaxy position as a reference point     
### NOTE: this is the same dir path that was provided to "mkcutouts.py"
#partone $fil1 $dirfitsname 
#
### Generates final color png images with an alpha layer of mask added to
### indicate the location of the simulated lensed images and compressed using "crush"
### algorithm
#parttwo $outdirname
#
### Copying all the generated FITS and PNGs in the output dir. 
### Note: local copies of the output files need to be manually deleted
#partthree $outdirname
#
