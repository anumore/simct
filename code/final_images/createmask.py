#!/usr/bin/env python
import numpy
import Image
import sys

if(len(sys.argv)!=4):
    print "./createmask.py finp.png mask.png out.png";
    sys.exit(0);
else:
    fname1=sys.argv[1];
    fname2=sys.argv[2];
    fout=sys.argv[3];

  
# Open image and mask using PIL
im = Image.open(fname1)
im_mask = Image.open(fname2)

# Convert to grayscale
im_mask = im_mask.convert('L')

# Get numpy array representation
mask = numpy.asarray(im_mask)

sz=len(mask[0]);

# Make array writable
mask.flags.writeable = True

# For SW purposes only need binary mask
mask[numpy.where(mask > 0)] = 255
 
win=5;

masknew=mask*1;

for ii in range(len(mask)):
    for jj in range(sz):
	if(mask[ii][jj]==255):
	    for ll in range(ii-win,ii+win+1):
		for kk in range(jj-win,jj+win+1):
		    if(ll<0 or ll>=len(mask) or kk<0 or kk>=sz):
			continue;
		    masknew[ll][kk]=255;

masknew[numpy.where(masknew != 255)] = 254
##im_mask = Image.fromarray(masknew)
##im_mask.save(fout)


# Convert to Image object and save
im_mask = Image.fromarray(masknew)

# Add mask to alpha of original image
im.putalpha(im_mask)

im.save(fout)

