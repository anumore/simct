#!/usr/bin/env python
from math import *
from scipy.interpolate import interp1d
from scipy.optimize import brentq
import numpy as np
import srclensprop as slp
import sigma_vd as svd
import scipy.integrate as sci
from input_qg import *

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
def setlogfile(lgfile,goutfile,qoutfile):
    global fp1,fp2,fp3;
    fp1=open(lgfile,"w");
    fp2=open(qoutfile,"w");
    fp3=open(goutfile,"w");


init_Phi_qso=0;
Phi_qso_spl=0;

init_P_gal=0;
P_gal_spl=0;

init_nsmlim_gal=0;
nsmlim_gal_spl=0;

########
## PART 1- For quasars
########

## The units here are h70^3 Mpc^{-3}/mag
def dPhibydM_qso(zz,mag):
    alp=-0.5; 
    kcorr=-2.5*(1 + alp)*log10(1+zz);
    Dlum=cc.Dlofz(zz)/p.hval;
    DM=5*log10(Dlum);
    Mabs=mag-kcorr-DM-25.0;
    zeta=2.98;
    xi=4.05;
    zstar=1.60;
  # hh=cc.H0/100.;
    if(zz<3):
	alpha=-3.31;
    else:
	alpha=-2.58;
    beta=-1.45;
    phistar=5.34E-6*pow(p.hval,3);
    fofz=(exp(zeta*zz)*(1+exp(xi*zstar)))/pow(np.sqrt(exp(xi*zz))+np.sqrt(exp(xi*zstar)),2);
    Mstar=-20.90+5*log10(p.hval)-2.5*log10(fofz);
    return phistar/(pow(10,(0.4*(alpha+1)*(Mabs-Mstar)))+pow(10,(0.4*(beta+1)*(Mabs-Mstar))));

## Units of Phi(z) will be h70^3 Mpc^{-3}
def initPhiqso():
  global Phi_qso_spl,init_Phi_qso;
  fp1.write("Initializing Phi(z)\n");
  zi=np.arange(zmin_qso,zmax_qso+0.01,0.01);
  Phii=zi*0.0;
  for ii in range(zi.size):
    Phii[ii]=sci.quad(lambda mag: dPhibydM_qso(zi[ii],mag),mmin_qso,mmax_qso)[0];
  Phi_qso_spl=interp1d(zi,Phii,kind='cubic');
  init_Phi_qso=1;

def dNbydz_qso(zz,vdisp,zred,q):
  ## \int dz dV/dz \int dM dPhi/dM(z) \sigma(zl,zs,q)
  global Phi_qso_spl,init_Phi_qso;
  if(init_Phi_qso==0):
    initPhiqso();
  ## Get Phi
  Phi=Phi_qso_spl(zz);

  ## Get cross-section in steradian
  bsis=slp.getreinst(zred,zz,vdisp);
  csect=slp.getcrosssect_num(bsis,q);

  ## Multiply by volume factor and return
  return Phi*1./cc.Eofz(zz)*(cspd/(p.hval*100))*(cc.Daofz(zz)/p.hval)**2.*csect*(1+zz)**2;

def findzqso(ztry,vdisp,zred,q,Ntarget):
  return Ntarget-sci.quad(lambda zz: dNbydz_qso(zz,vdisp,zred,q),zmin_qso,ztry)[0];

def findmagqso(magtry,zsrc,Phitarget):
  result=Phitarget-sci.quad(lambda mag: dPhibydM_qso(zsrc,mag),mmin_qso,magtry)[0];
  return result;

def Nsrc_qso(magg,magr,zred,q):
  global boost_csect_qso
  ## Get the velocity dispersion
  vdisp=svd.getsigma(magg,magr,zred);
  
  bsist=slp.getreinst(zred,zmax_qso,vdisp)*206264.8;

  if (bsist<1.2):
      return 0,0,0,0,0.;
  else:
      ## Now calculate the number of quasars around this object within its
      ## lensing cross-section. 
      ## This is given by the following integral
      ## \int dz dV/dz \int dM dPhi/dM(z) \sigma(zl,zs,q)
      ##
      ## For quasars, let us integrate from zmin=1.0 to 5.0
      Nsrcmean=sci.quad(lambda zz: dNbydz_qso(zz,vdisp,zred,q),zmin_qso,zmax_qso)[0];
      Nsrcmean_boost=Nsrcmean*boost_csect_qso;

      ## Return a Poisson deviate
      Nreal=np.random.poisson(Nsrcmean_boost);

      fp1.write("This lens has %f quasars behind it on average and Poisson deviate is %d\n"%(Nsrcmean_boost,Nreal));

      listmag=[];
      listz=[];
      rands=[];

      for ii in range(Nreal):
	  qsozsrc=0.0;
	  while(qsozsrc<zred):
	      rr=np.random.random();
	      Ntarg=rr*Nsrcmean;
	      ## Need a root finder to get the redshift
	      qsozsrc=brentq(findzqso,zmin_qso,zmax_qso,args=(vdisp,zred,q,Ntarg),xtol=1.e-3);
	  rands=np.append(rands,rr);
	  listz=np.append(listz,qsozsrc);
	  done=False;
	  while(not done):
	      try:
		  rr=np.random.random();
		  Phitarg=rr*Phi_qso_spl(qsozsrc);
		  rands=np.append(rands,rr);
		  qsomag=brentq(findmagqso,mmin_qso,mmax_qso,args=(qsozsrc,Phitarg),xtol=1.e-3);
		  done=True;
	      except ValueError:
		  print "checking", findmagqso(mmin_qso,qsozsrc,Phitarg), findmagqso(mmax_qso,qsozsrc,Phitarg),mmin_qso, mmax_qso, Phitarg;
		  print "Done is:",done;
	  listmag=np.append(listmag,qsomag);

      return Nreal,listmag,listz,rands,vdisp;

## Calculate the total number of qsos in a survey
def Nqso():
  global Phi_qso_spl,init_Phi_qso;
  if(init_Phi_qso==0):
    initPhiqso();
  
  def dNqso(zz):
      ## Get Phi
      Phi=Phi_qso_spl(zz);
      return Phi*1./cc.Eofz(zz)*(cspd/(p.hval*100))*(cc.Daofz(zz)/p.hval)**2.*(1+zz)**2;

  ## Units are steradian^(-1) 
  totNqso=sci.quad(lambda zz:dNqso(zz),zmin_qso,zmax_qso)[0];


########
## PART 2- For galaxies
########

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


def dNbydz_gal(zz,vdisp,zred,q):
  ## \int dz dV/dz p(zs,mlim) \sigma(zl,zs,q)
  global P_gal_spl,init_P_gal;
  if(init_P_gal==0):
    initPhigal(mmax_gal);
  ## Get Phi
  Phi=P_gal_spl(zz);

  ## Get cross-section in steradian
  bsis=slp.getreinst(zred,zz,vdisp);
  csect=slp.getcrosssect_num(bsis,q);

  return Phi*csect;

def findzgal(ztry,vdisp,zred,q,Ntarget):
  return Ntarget-sci.quad(lambda zz: dNbydz_gal(zz,vdisp,zred,q),zmin_gal,ztry)[0];

def findmaggal(magtry,zsrc,Phitarget):
  z0=0.13*magtry-2.2;
  beta=1.5;
  Pmagzs=beta*(zsrc/z0)**2.0*exp(-(zsrc/z0)**beta)/z0;
  return Phitarget-nsmlim_gal_spl(magtry)*Pmagzs;

def Nsrc_gal(magg,magr,zred,q):
  global boost_csect_gal,init_nsmlim_gal;

  listmag=[];
  listz=[];
  rands=[];

  if(zred>zmax_gal):
      return 0,listmag,listz,rands;

  ## Initialize ns(mlim) if not done before
  if(init_nsmlim_gal==0):
      initnsmlim();

  ## Get the velocity dispersion
  vdisp=svd.getsigma(magg,magr,zred);
  bsist=slp.getreinst(zred,zmax_gal,vdisp)*206264.8;
  
  if (bsist<1.2):
      return 0,0,0,0,0.;
  
  else:
      ## Calculate the number of galaxies around this object within its
      ## lensing cross-section. 
      ## This is given by the following integral
      ## \int dz dV/dz P(zs,mlim) \sigma(zl,zs,q)
      #
      Nsrcmean=sci.quad(lambda zz: dNbydz_gal(zz,vdisp,zred,q),zmin_gal,zmax_gal)[0];
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
	      galzsrc=brentq(findzgal,zmin_gal,zmax_gal,args=(vdisp,zred,q,Ntarg),xtol=1.e-3);
	  rands=np.append(rands,rr);
	  listz=np.append(listz,galzsrc);
	  done=False;
	  while(not done):
	      try:
		  rr=np.random.random();
		  Phitarg=rr*P_gal_spl(galzsrc)*nsmlim_gal_spl(mmax_gal);
		  rands=np.append(rands,rr);
		  galmag=brentq(findmaggal,mmin_gal,mmax_gal,args=(galzsrc,Phitarg),xtol=1.e-3);
		  done=True;
	      except ValueError:
		  print "checking", findmaggal(mmin_gal,galzsrc,Phitarg), findmaggal(mmax_gal,galzsrc,Phitarg),mmin_gal, mmax_gal, Phitarg;
		  print "Done is:",done;
	  listmag=np.append(listmag,galmag);
	  

      return Nreal,listmag,listz,rands,vdisp;

