#!/usr/bin/env python
from math import *;
from scipy.interpolate import interp1d;
from scipy.optimize import brentq;
import numpy as np;
import genimg as gn
import scipy.integrate as sci
from input_g import *

###########################################################
## PURPOSE:
## Use source redshift distributions and luminosity functions to calculate the
## source properties. Calculate the source number density and lens cross-section
## to derive the number of sources per lens along with the source redshift and
## i-band magnitude 
###########################################################

fp1=0;
fp2=0;
fp3=0;
def setlogfile(lgfile,goutfile,gcfile):
    global fp1,fp2,fp3;
    fp1=open(lgfile,"w");
    fp2=open(goutfile,"w");
    fp3=open(gcfile,"w");

## Input: Magnitude of galaxy in the g and r band in Megacam filters and its
## redshift and its axis ratio
##
## Output: The number of galaxies this cluster can lens

init_P_gal=0;
P_gal_spl=0;
init_nsmlim_gal=0;
nsmlim_gal_spl=0;


## This gives a dimensionless P(zs,mlim)
def initPhigal(mlim):
  global P_gal_spl,init_P_gal;
  fp1.write("Initializing Phi(z)\n");
  zi=np.arange(zmin_gal,zmax_gal+0.01,0.01);
  Phii=np.arange(zi.size)*0.0;
  z0=0.13*mlim-2.2;
  beta=1.5;
  for ii in range(zi.size):
      Phii[ii]=beta*(zi[ii]/z0)**2.0*exp(-(zi[ii]/z0)**beta)/z0;
  P_gal_spl=interp1d(zi,Phii,kind='cubic');
  init_P_gal=1;

def dnsbydm(mag):
  ## deg^-2 --> radian^-2
  n0=3e3*(180./pi)**2.0;
  m1=20.;
  a=0.30;
  b=0.56;
  return n0/sqrt(10.0**(2*a*(m1-mag))+10.**(2*b*(m1-mag)) )

def initnsmlim():
  global nsmlim_gal_spl,init_nsmlim_gal;
  fp1.write("Initializing ns(mlim)\n");
  magi=np.arange(mmin_gal,mmax_gal+0.01,0.01);
  nsi=np.arange(magi.size)*0.;
  for ii in range(magi.size):
      nsi[ii]=sci.quad(lambda xmag:dnsbydm(xmag),mmin_gal,magi[ii])[0];
  nsmlim_gal_spl=interp1d(magi,nsi,kind='cubic');
  init_nsmlim_gal=1;


def dNbydz_gal(zz,kappa1,rscale,vdisp,gpx,gpy,ell,ellpa,zred,cnt,num):
  ## \int dz dV/dz p(zs,mlim) \sigma(zl,zs,q)
  global P_gal_spl,init_P_gal;
  if(init_P_gal==0):
    initPhigal(mmax_gal);
  ## Get Phi
  Phi=P_gal_spl(zz);
  
  if(gn.init_cs_spl==0):
    gn.init_crosssect(kappa1,rscale,vdisp,gpx,gpy,ell,ellpa,zred,cnt,pixsc,num,zmin_gal,zmax_gal);
  
  ## Then get cross-section in steradian
  csect=gn.cross_sect_spl(zz);
  return Phi*csect;

def findzgal(ztry,kappa1,rscale,vdisp,gpx,gpy,ell,ellpa,zred,cnt,num,Ntarget):
  return Ntarget-sci.quad(lambda zz: dNbydz_gal(zz,kappa1,rscale,vdisp,gpx,gpy,ell,ellpa,zred,cnt,num),zmin_gal,ztry)[0];

def findmaggal(magtry,zsrc,Phitarget):
  z0=0.13*magtry-2.2;
  beta=1.5;
  Pmagzs=beta*(zsrc/z0)**2.0*exp(-(zsrc/z0)**beta)/z0;
  return Phitarget-nsmlim_gal_spl(magtry)*Pmagzs;

def Nsrc_gal(kappa1,rscale,vdisp,gpx,gpy,ell,ellpa,zred,cnt,num):
  global boost_csect_gal,init_nsmlim_gal;

  listmag=[];
  listz=[];
  rands=[];

  if(zred>zmax_gal):
      return 0,listmag,listz,rands;

  ## Initialize ns(mlim) if not done before
  if(init_nsmlim_gal==0):
      initnsmlim();
  
  ## Now calculate the number of galaxies around this object within its
  ## lensing cross-section. 
  ## This is given by the following integral
  ## \int dz dV/dz P(zs,mlim) \sigma(zl,zs,q)
  ##
  ## For galaxies, let us integrate from zmin=1.0 to 3.0
  
  gn.init_cs_spl=0;
  Nsrcmean=sci.quad(lambda zz: dNbydz_gal(zz,kappa1,rscale,vdisp,gpx,gpy,ell,ellpa,zred,cnt,num),zmin_gal,zmax_gal)[0];
  Nsrcmean_boost=Nsrcmean*nsmlim_gal_spl(mmax_gal)*boost_csect_gal;

  ## Return a Poisson deviate
  Nreal=np.random.poisson(Nsrcmean_boost);

  fp1.write("This lens has %f galaxies behind it on average and Poisson deviate is %d\n"%(Nsrcmean_boost,Nreal));
  
  for ii in range(Nreal):
      galzsrc=0.0;
      trials=0
      while(galzsrc<zred):
          rr=np.random.random();
          Ntarg=rr*Nsrcmean;
          ## Need a root finder to get the redshift
          galzsrc=brentq(findzgal,zmin_gal,zmax_gal,args=(kappa1,rscale,vdisp,gpx,gpy,ell,ellpa,zred,cnt,num,Ntarg));
      rands=np.append(rands,rr);
      listz=np.append(listz,galzsrc);
      done=False;
      while(not done):
          try:
	      rr=np.random.random();
	      Phitarg=rr*P_gal_spl(galzsrc)*nsmlim_gal_spl(mmax_gal);
	      rands=np.append(rands,rr);
	      galmag=brentq(findmaggal,mmin_gal,mmax_gal,args=(galzsrc,Phitarg));
	      done=True;
          except ValueError:
	      print "checking", findmaggal(mmin_gal,galzsrc,Phitarg), findmaggal(mmax_gal,galzsrc,Phitarg),mmin_gal, mmax_gal, Phitarg;
	      print "Done is:",done;
      listmag=np.append(listmag,galmag);

  return Nreal,listmag,listz,rands;

