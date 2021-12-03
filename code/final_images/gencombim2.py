#!/usr/bin/env python
# reads ra,dec, field name and id
# calculates in which cutout the candidate is located and writes the file for making color png
#===============

import numpy as np;
from math import *;
from subprocess import call;
import sys


if(len(sys.argv)!=3):
    print "./gencombim2.py dirname inpcat";
    sys.exit(0);
else:
    dirname=sys.argv[0];
    inpcat=sys.argv[1];

## Input Catalog
gid,gra,gdec,fld,fldid,gxpx,gypx=np.loadtxt(inpcat,dtype={'names': ('gid','gra', 'gdec','fld','fldid','gxpx','gypx'), 'formats': ('d', 'f', 'f', 'S13', 'S3','f','f')},unpack=True);
## List with xmin,xmax,ymin,ymax of the contiguous cutouts made for each survey
## tile
offsetx,offsety=np.loadtxt("cutoutlist",usecols=(0,2),unpack=True);
################

band=['u','g','r','i','z'];
Ntot=len(band);

call("rm imdir/lstfile* imdir/rm*_* imdir/rimst_*",shell=1);

def writefiles(dirname,var,var1,fldid,gid,ofstx,ofsty):
    for kk in range(Ntot):
	fp1.write('cp %s/%s/%s/%s_%s.fits ../outfits/%s_%s_t.fits \n'%(dirname,fldid,band[kk],var,band[kk],var,band[kk]));
	call("awk \'{ printf(\"../gout/imoutp_'%s'_'%d'_'%s'.fits,'%f','%f' \\n \" );}\' gg >> imdir/lstfile_'%s'_'%s'"%(fldid,gid,band[kk],ofstx,ofsty,var1,band[kk]),shell=1);
	fp.write('fimgmerge ../outfits/%s_%s_t.fits  @lstfile_%s_%s ../outfits1/%s_%s.fits \n'%(var,band[kk],var1,band[kk],var,band[kk]));
	fp2.write('fimgmerge blnkcutout.fits  @lstfile_%s_%s ../blanksims/%s_%s.fits \n'%(var1,band[kk],var,band[kk]));
	fp3.write('./getimstat.py ../outfits1/%s_%s.fits\n'%(var,band[kk]));


fpw=open('imdir/rmrgall','w');
fpw1=open('imdir/rcpyall','w');
fpw2=open('imdir/rmrgball','w');
fpw3=open('imdir/rimsall','w');
fpw4=open('imdir/idpxlst','w');

fpw.write('#!/bin/bash \n');
fpw1.write('#!/bin/bash \n');
fpw2.write('#!/bin/bash \n');
fpw3.write('#!/bin/bash \n');

flag=1;
for ii in range(gid.size):
    ix=int((gxpx[ii]+0.5)/386.0);
    iy=int((gypx[ii]+0.5)/386.0);
    cutoutno=(ix-1)*50+iy;
    modx=(gxpx[ii])%386.+1.;
    mody=(gypx[ii])%386.+1.;
    dgx=(gxpx[ii]-101);
    dgy=(gypx[ii]-101);
    ptx1=-9e5;
    ptx2=-9e5;
    pty1=-9e5;
    pty2=-9e5;
    
    if (modx>54):
       ptx1=(ix)*50;
       ptx2=-9e5;
    
    if (mody>54):
       pty1=iy+1
       pty2=-9e5;
    
    if (ix>0 and (modx<=54 and modx>=1.) and ix<50):
       ptx1=(ix-1)*50;
       ptx2=ix*50;
    
    if (iy>0 and (mody<=54 and mody>=1.) and iy<50):
       pty1=iy
       pty2=iy+1;

    if (ix==50 and (modx<=54 and modx>=1.)):
       ptx1=(ix-1)*50;
       ptx2=-9e5;
    
    if (iy==50 and (mody<=54 and mody>=1.)):
       pty1=iy
       pty2=-9e5;
    
    
    ctn1=ptx1+pty1;
    ctn2=ptx1+pty2;
    ctn3=ptx2+pty1;
    ctn4=ptx2+pty2;
    if(flag==1):
        fp=open('imdir/rmrg_%s'%(fldid[ii]),'w');
        fp1=open('imdir/rmgcpy_%s'%(fldid[ii]),'w');
        fp2=open('imdir/rmrgb_%s'%(fldid[ii]),'w');
        fp3=open('imdir/rimst_%s'%(fldid[ii]),'w');
        
        fp.write('#!/bin/bash \n');
        fp1.write('#!/bin/bash \n');
        fp2.write('#!/bin/bash \n');
        fp3.write('#!/bin/bash \n');
        
        fpw.write('./rmrg_%s \n'%(fldid[ii]));
        fpw1.write('./rmgcpy_%s \n'%(fldid[ii]));
        fpw2.write('./rmrgb_%s \n'%(fldid[ii]));
        fpw3.write('./rimst_%s \n'%(fldid[ii]));
        
        flag=0;
    
    
    if(fldid[ii]==fldid[ii-1] or flag==0):
	print "Writing:",ii+1,int(gid[ii]),fldid[ii];
	if(ctn1>0):
	   var="CFHTLS_%s_%04d"%(fldid[ii],ctn1);
	   var1="%s_%04d"%(fldid[ii],ctn1);
	   ofstx=dgx-offsetx[ctn1-1];
	   ofsty=dgy-offsety[ctn1-1];
	   fpw4.write(' %s %d %f %f \n'%(var1,gid[ii],ofstx+102,ofsty+102));
	   writefiles(dirname,var,var1,fldid[ii],gid[ii],ofstx,ofsty);
	
	if(ctn2>0):
	   var="CFHTLS_%s_%04d"%(fldid[ii],ctn2);
	   var1="%s_%04d"%(fldid[ii],ctn2);
	   ofstx=dgx-offsetx[ctn2-1];
	   ofsty=dgy-offsety[ctn2-1];
	   fpw4.write(' %s %d %f %f \n'%(var1,gid[ii],ofstx+102,ofsty+102));
	   writefiles(dirname,var,var1,fldid[ii],gid[ii],ofstx,ofsty);
	
	if(ctn3>0):
	   var="CFHTLS_%s_%04d"%(fldid[ii],ctn3);
	   var1="%s_%04d"%(fldid[ii],ctn3);
	   ofstx=dgx-offsetx[ctn3-1];
	   ofsty=dgy-offsety[ctn3-1];
	   fpw4.write(' %s %d %f %f \n'%(var1,gid[ii],ofstx+102,ofsty+102));
	   writefiles(dirname,var,var1,fldid[ii],gid[ii],ofstx,ofsty);
	
	if(ctn4>0):
	   var="CFHTLS_%s_%04d"%(fldid[ii],ctn4);
	   var1="%s_%04d"%(fldid[ii],ctn4);
	   ofstx=dgx-offsetx[ctn4-1];
	   ofsty=dgy-offsety[ctn4-1];
	   fpw4.write(' %s %d %f %f \n'%(var1,gid[ii],ofstx+102,ofsty+102));
	   writefiles(dirname,var,var1,fldid[ii],gid[ii],ofstx,ofsty);
	fp.write("echo \"Done: %d %d \" \n"%(ii+1,gid[ii]));
	fp1.write("echo \"Done: %d %d \" \n"%(ii+1,gid[ii]));
	fp2.write("echo \"Done: %d %d \" \n"%(ii+1,gid[ii]));
	fp3.write("echo \"Done: %d %d  \"\n"%(ii+1,gid[ii]));

    if (ii==gid.size-1 or fldid[ii]!=fldid[ii+1] ):	    
        fp.write("echo \"Done:rmrg_%s \" \n"%(fldid[ii]));
        fp1.write("echo \"Done:rmrgcpy_%s \" \n"%(fldid[ii]));
        fp2.write("echo \"Done:rmrgb_%s \" \n"%(fldid[ii]));
        fp3.write("echo \"Done:rimst_%s  \"\n"%(fldid[ii]));
	fp.close();
        fp1.close();
        fp2.close();
        fp3.close();
        flag=1;

fpw.close();	
fpw1.close();	
fpw2.close();	
fpw3.close();	
fpw4.close();	
call("chmod +x imdir/rmrg_* imdir/rmrgb_* imdir/rmgcpy_* imdir/rimst_* imdir/rmrgall imdir/rcpyall imdir/rmrgball imdir/rimsall",shell=1);
