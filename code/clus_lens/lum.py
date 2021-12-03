#!/usr/bin/env python
import sys
import numpy as np
from math import *
from scipy.interpolate import interp1d
from scipy.integrate import quad
from input_g import *

###########################################################
## INPUTS: 
## FileFor_Kapp_s_R_s.dat which has Redshift, M200, Mvir, Rs, Rhos (see below) 
## 
## PURPOSE:
## Calculate the rho_s and r_s values for a halo given the BCG magnitude and
## redshift
###########################################################


## Values from Behroozi et al. 2012
eps0=-1.785; 
epsa=-0.496;
epsz=-0.167;
epsa2=-0.101;

xM10=11.516;
xM1a=-1.894;
xM1z=-0.368;

alp0=1.432;
alpa=-0.088;

delta0=3.581;
deltaa=4.567;
deltaz=1.147;

gamma0=0.349;
gammaa=1.484;
gammaz=0.243;

readfile=0;

zred=0;
m200=0;
mvir=0;
rs=0;
rhos=0;

def readfileforksrs():

## Columns in FileFor_Kapp_s_R_s.dat
## redshift
## log10(M200): in hinv Msun
## log10(Mvir): in hinv Msun
## log10(Rs): in hinv Mpc
## log10(Rhos): in hinv Msun/(hinv Mpc)^3
    global zred,m200,mvir,rs,rhos,readfile;
    zred,m200,mvir,rs,rhos=np.loadtxt(fileforkappa,unpack=1);
    readfile=1;

def gammaofz(zz):
    global alp0,alpa,delta0,deltaa,deltaz,gamma0,gammaa,gammaz;
    a=1./(1.+zz);
    nu=exp(-4*a*a);
    gamma=gamma0+(gammaa*(a-1.)+gammaz*zz)*nu;
    return gamma;

def fbeh12(x,zz):
    global alp0,alpa,delta0,deltaa,deltaz,gamma0,gammaa,gammaz;
    a=1./(1.+zz);
    nu=exp(-4*a*a);
    alp=alp0+alpa*(a-1.0)*nu;
    delta=delta0+(deltaa*(a-1.)+deltaz*zz)*nu;
    gamma=gamma0+(gammaa*(a-1.)+gammaz*zz)*nu;
    return -np.log10(10.**(-alp*x)+1)+delta*(np.log10(1.+np.exp(x)))**gamma/( 1. + np.exp(10.**-x) )

def SHMRbeh(xmh,zz):
    global eps0,epsa,epsa2,epsz,xM10,xM1a,xM1z;
    a=1./(1.+zz);
    nu=exp(-4*a*a);
    eps=eps0+(epsa*(a-1.0)+epsz*zz)*nu+epsa2*(a-1.) 
    xM1=xM10+(xM1a*(a-1.0)+xM1z*zz)*nu;
    xmstel=eps+xM1+fbeh12(xmh-xM1,zz)-fbeh12(0.0,zz);
    return xmstel;

def pmsnmdm(xmstel,x,zz,getmstel):
    m=10.0**x;
    fac=cc.nofm(m,zz)*m*log(10.);
    sig=0.18;
    return exp(-(xmstel-getmstel(x))**2.0/2./sig**2)*fac;

def fetchMh(xmstel,getmstel, zz):
    ## First find the mean and the sigma of P(log M|log M*)
    num,err=quad(lambda x: pmsnmdm(xmstel,x,zz,getmstel)*x ,9.5,15.9);
    denom,err=quad(lambda x: pmsnmdm(xmstel,x,zz,getmstel) ,9.5,15.9);
    num2,err=quad(lambda x: pmsnmdm(xmstel,x,zz,getmstel)*x**2 ,9.5,15.9);
    sigmalogm=sqrt(num2/denom-(num/denom)**2.);
    avlogm=num/denom;
    #return avlogm,sigmalogm,np.random.normal(avlogm,sigmalogm);
    return np.random.normal(avlogm,sigmalogm);
    
def getRhos_Rs_i(mi,zz):
    global readfile,cc;
    ## in Mpc
    #Dlum=cc.DlumMpc(zz);
    
    ## hval=0.72;
    ## in hinvMpc
    Dlum=cc.Dlofz(zz);

    ## Dlum should be in hinv Mpc
    mi=mi-5.0*log10(Dlum/p.hval)-25.0

    ## Get the luminosity in terms of Lsun
    xlum=-0.4*(mi-4.54);

    ## Convert to M_*/Msun assuming a constant mass to light ratio of 3*0.7 h Msun/Lsun
    xmstel=xlum+log10(2.1);
    sys.stdout.flush();

    ## Get the mstel_beh-mhalo_beh relation at the redshift
    xMh_beh=np.arange(9.0,16.0,0.1);
    xmst_beh=SHMRbeh(xMh_beh,zz);

    getmstel=interp1d(xMh_beh,xmst_beh);

    ## Adds in a random scatter
    xMhalo=fetchMh(xmstel,getmstel,zz);

    ## Get Mhalo in hinv Msun
    xMhalo=xMhalo+log(p.hval);

    ## Get kappas and rs
    if(readfile==0):
      readfileforksrs();

    ## Find a redshift and halo mass match
    idx=np.argmin(np.fabs(zred-zz)+np.fabs(xMhalo-mvir));
    rsmatch=10.0**rs[idx];
    rhosmatch=10.0**rhos[idx];

    return rsmatch,rhosmatch

if __name__ == '__main__':
    xMh_beh=np.arange(9.0,16.0,0.1);
    zz=0.5
    xmst_beh=SHMRbeh(xMh_beh,zz);
    getmstel=interp1d(xMh_beh,xmst_beh);

    xmstel=np.arange(10.0,12.0,0.1);
    for i in range(xmstel.size):
	xMhalo=fetchMh(xmstel[i],getmstel,zz);
	##avMh,sigMh,xMhalo=fetchMh(xmstel[i],getmstel,zz);
	#print xmstel[i],xMhalo;
	#print xmstel[i],avMh;
