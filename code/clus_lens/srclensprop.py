#!/usr/bin/env python
from math import *
from pylab import *
from scipy.interpolate import interp1d
from scipy.integrate import quad
import numpy as np
from input_g import *

###########################################################
## PURPOSE:
## This module will provide a random source position which lies within
## the caustic of a singular isothermal ellipsoid. The definition
## conventions are based upon Keeton Mao and Witt 2004. I set the
## external shear to be equal to zero in their formulae.
##
## INPUTS:
## b_I: This is the einstein radius of the SIE
##      b_I = b_SIS eps_3/sin^{-1}(eps_3), where eps_3 is the
##      eccentricity of the mass distribution
##
## q: Projected axis ratio
##      q=sqrt(q3^2sin^2i+cos^2i) where q3 is the 3d acis ratio of the
##      ellipsoid and i is inclination angle
##
##      eps_3=(1-q3^2)^{1/2}
##
## pa: Position angle in degrees
##
##
## Output:
## srcx,srcy : The x and y positions of the sources where the major axis is
## aligned with the x axis.
##
## Note that we also provide a wrapper (srcpos_bsisq) which uses a value for b_SIS
## and q to define the projected surface mass density. Comparing Keeton, Mao and
## Witt 2000 equations with Kormann et al. 94 yields:
##
## b_I = sqrt(q)*b_{SIS}
##
## Although it seems that magically the inclination angle dependence vanishes
## out from kappa. This equation is on average ok given the oblateness and
## prolateness of the halos (Oguri, priv comm., see also paper by Chae 2003).
## 
###########################################################

## Constants
fid_b_I=10.0;
init_cs_spl=0;
cross_sect_spl=0;

def getcaustics(b_I,q):
  ## Do not change the next line or face consequences
  t=np.arange(0.0,pi/2.0+pi/100.0,pi/200.0);

  den=(1+q**2)-(1-q**2)*np.cos(2*t)  
  r=sqrt(2.)*b_I/np.sqrt(den)
  x=r*np.cos(t);
  y=r*np.sin(t);

  ## Generate parametric function for the tangential caustic
  xi=np.sqrt(2.*(1-q**2)/den);
  u=np.cos(t)*r;
  v=np.sin(t)*r;
  if(q!=1.0):
      u=u-b_I*np.arctan(xi*np.cos(t))/sqrt(1-q**2);
      v=v-b_I*np.arctanh(xi*np.sin(t))/sqrt(1-q**2);

  ## Generate parametric function for the radial caustic
  up=b_I*0.0;
  vp=b_I*0.0;
  if(q!=1.0):
      up=-b_I*np.arctan(xi*np.cos(t))/sqrt(1-q**2);
      vp=-b_I*np.arctanh(xi*np.sin(t))/sqrt(1-q**2);

  ## Set up splines for the tangential and radial caustic
  r1=np.sqrt(u**2+v**2);
  r2=np.sqrt(up**2+vp**2);
  
  t1=[0.]*u.size; 
  t2=[0.]*up.size; 
  for ii in range(u.size):
      if(u[ii]==0.):
	  t1[ii]=np.pi/2.;
      else:
	  t1[ii]=arctan(np.abs(v[ii]/u[ii]));
      if(t[ii]>pi/2.0):
	  t1[ii]=pi-t1[ii];
  
  for ii in range(up.size):
      if(up[ii]==0.):
	  t2[ii]=np.pi/2.;
      else:
	  t2[ii]=arctan(np.abs(vp[ii]/up[ii]));
      if(t[ii]>pi/2.0):
	  t2[ii]=pi-t2[ii];
  c1=interp1d(t1,r1,kind='cubic');
  c2=interp1d(t2,r2,kind='cubic');

  ## Calculate maximum radial distance between the caustics and the origin
  maxr=np.max([r1,r2]);

  return c1,c2,maxr;

## Given the sie caustics, get the source position
def srcpos_sie(b_I,q,pa):
  ## Generate parametric theta values and the caustics at that theta
  c1,c2,maxr=getcaustics(b_I,q);

  ## Generate random points within a circle with radius maxr until you find one
  ## inside the caustic
  while(1):
    rgen=sqrt(np.random.rand(1))*maxr;
    tgen=np.random.rand(1)*(pi/2.);

    if(rgen<c1(tgen) or rgen<c2(tgen) ):
      break;

  ## Generate the positions in the first quadrant
  xret,yret=rgen*cos(tgen),rgen*sin(tgen);

  ## Assign quadrant
  quadr=np.random.rand(1);
  if(quadr<0.5):
    xret=-xret;
  quadr=np.random.rand(1);
  if(quadr<0.5):
    yret=-yret;

  ## Rotate by pa
  rad=pa*pi/180.0;
  xnew=xret*cos(rad)-yret*sin(rad)
  ynew=xret*sin(rad)+yret*cos(rad)

  return xnew,ynew;

def srcpos_bsisq(bsis,q,pa):
  xret,yret=srcpos_sie(bsis*sqrt(q),q,pa);
  return xret,yret;

## Gives cross-section in units of bsis
def getcrosssect(bsis,q):
  b_I=bsis*sqrt(q);

  ## Generate parametric theta values and the critical curve
  t=np.arange(0.0,pi/2.0+pi/200.0,pi/200.0);
  c1,c2,maxr=getcaustics(b_I,q);
  r1=c1(t);
  r2=c2(t);

  for i in range(r1.size):
    r1[i]=max(r1[i],r2[i]);

  ## We want to integrate \int_0^{2pi} 1/2 r^2(t) dt= 4 \int_0^{pi/2} 1/2 r^2(t) dt
  cnew=interp1d(t,4*0.5*r1*r1);
  
  #Now integrate under the curve
  result,err=quad(lambda x:cnew(x),0.,pi/2.,epsrel=1.0E-3)
  ##print result;

  #This is the angular area is in units^2 where units is the unit of b. 
  return result;
   
## Einstein radius in radians
def getreinst(zlens,zsrc,sigma):
    Ds=cc.Daofz(zsrc);
    Dds=cc.Daofzlh(zlens,zsrc);
    reinst=4*pi*(sigma/cspd)**2.*(Dds/Ds);
    return reinst;

def getreinst2(zlens,zsrc,sigma):
    Ds=cc.Daofz(zsrc);
    Dds=cc.Daofzlh(zlens,zsrc);
    reinst=4*pi*(sigma/cspd)**2.*(Dds/Ds);
    Drat=Ds/Dds;
    return reinst,Drat;

## Module to plot critical curves and caustics
def plotcaust(b_I,q,pa):
    t=np.arange(0.0,2*pi,0.0001);
    den=(1+q**2)-(1-q**2)*np.cos(2*t)  
    r=sqrt(2.)*b_I/np.sqrt(den)
    x=r*np.cos(t);
    y=r*np.sin(t);
    rad=pa*pi/180.;
    xnew=x*cos(rad)-y*sin(rad)
    ynew=x*sin(rad)+y*cos(rad)
    plot(xnew,ynew);
    
    xi=np.sqrt(2.*(1-q**2)/den);
    u=np.cos(t)*r-b_I*np.arctan(xi*np.cos(t))/sqrt(1-q**2);
    v=np.sin(t)*r-b_I*np.arctanh(xi*np.sin(t))/sqrt(1-q**2);
    
    up=-b_I*np.arctan(xi*np.cos(t))/sqrt(1-q**2);
    vp=-b_I*np.arctanh(xi*np.sin(t))/sqrt(1-q**2);
    
    unew=u*cos(rad)-v*sin(rad)
    vnew=u*sin(rad)+v*cos(rad)
    plot(unew,vnew);
    upnew=up*cos(rad)-vp*sin(rad)
    vpnew=up*sin(rad)+vp*cos(rad)
    plot(upnew,vpnew);

## First note that the cross-section in sie just depends upon q where as the
## dependence on b_I can be easily predicted. As source redshift changes, only
## b_I changes. This implies we can initialize a cross-section spline which is a
## function of only q for a fiducial b_I.

def init_crosssect():
  global cross_sect_spl,init_cs_spl;
  print "Initializing Cross-section"
  ## Use the fiducial b_I, crosssect requires me to pass bsis.
  if(0):
      qi=np.arange(0.1,1.0,0.001);
      csecti=qi*0.0;
      for ii in range(qi.size):
	print qi[ii];
	csecti[ii]=getcrosssect(fid_b_I/sqrt(qi[ii]),qi[ii]);
      np.savetxt("Crosssect.dat",np.transpose([qi,csecti]));
  else:
      qi,csecti=np.loadtxt("Crosssect.dat",unpack=1);
  cross_sect_spl=interp1d(qi,csecti,kind='cubic');
  init_cs_spl=1;

def getcrosssect_num(bsis,q):
  global cross_sect_spl,init_cs_spl;
  if(q>0.999):
      q=0.999;
  if(q<0.1):
      q=0.1;
  if(init_cs_spl==0):
    init_crosssect();
  b_I=bsis*sqrt(q);
  return cross_sect_spl(q)*(b_I/fid_b_I)**2;

## Example program which generates 1000 points within the caustics and plots them
#x=np.arange(0.0,1.0E3,1.0)*0;
#y=np.arange(0.0,1.0E3,1.0)*0;
#for i in range(x.size):
#  print i;
#  x[i],y[i]=srcpos_sie(2.0,0.2,30.0);
#
#plotcaust(2.0,0.2,30.0);
#scatter(x,y,s=0.5);
#show();


