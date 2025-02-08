#!/usr/bin/env python
from math import *
from subprocess import call
import numpy as np
import sys
import random
import os    
#from astropy.cosmology import FlatLambdaCDM
#cosmo = FlatLambdaCDM(H0=72,Om0=0.26)

#myseed = 72398
#np.random.seed(myseed)


## Extract source position and some lensed image properties for quasars
def srcposrnq(reinst,ell,ell_pa,sh_str,sh_pa,glno,flag,gid,lenscode):
    limr=reinst.size
    if(flag):
        fp=open('mckq/gllens%d.in'%(glno),'w')
 
 #      np.savetxt(fp,['set rscale=15'],fmt='%s');
 #      np.savetxt(fp,['set ngrid1=40'],fmt='%s');
 #      np.savetxt(fp,['set ngrid2=40'],fmt='%s');
 #      np.savetxt(fp,['set verbose=0'],fmt='%s');
 #      np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');
        fp.write("set rscale=15 \n set ngrid1=40 \n set ngrid2=40 \n set verbose=0 \n set omitcore=0.02 \n")
 
        for ii in range(limr):
         #  np.savetxt(fp,['startup 1 1'],fmt='%s');
            fp.write("startup 1 1\n")
            fp.write('alpha %f 0 0 %f %f %f %f 0 0 1 \n' %(reinst[ii],ell,ell_pa,sh_str,sh_pa))
            fp.write(" 0 0 0 0 0 0 0 0 0 0 \n")
         #  np.savetxt(fp,[' 0 0 0 0 0 0 0 0 0 0'],fmt='%s');
            fp.write("mock1 mckq/out%d_%d 500 2 \n \n" %(glno,ii+1))
 
        #np.savetxt(fp,['quit'],fmt='%s');
        fp.write("quit")
        fp.close()
     
    print("Running mckq/gllens%d.in"%(glno))
    call('%s mckq/gllens%d.in > mckq/gltmpout_%d'%(lenscode,glno,glno),shell=1)
    
    ii=0
    srcx=np.arange(limr)*0
    srcy=np.arange(limr)*0
    smag=np.arange(limr)*0
    summ=np.arange(limr)*0
    imno=np.arange(limr)*0
     
    for ii in range(limr):
        ip=ii+1
        #np.random.seed(myseed)
        imgnos=np.random.uniform(0,1)
        ## Loop to select a source that has 2 or more lensed images such that a
        ## double is chosen 50% of the times
        cond=1
        while(cond<3):
            if(imgnos<0.5):
                call("awk -f scr2img.awk mckq/out'%d'_'%d' > mckq/tp_%d"%(glno,ip,glno),shell=True);
                ## if no doubly images sources (unlikely to happen), try sources with no. of lensed
                ## images >2
                if( os.path.getsize("mckq/tp_%d"%(glno))<=0):
                    imgnos=0.6
                    cond=cond+1
                else:
                    rwno=np.loadtxt("mckq/tp_%d"%(glno),unpack=True)
                    cond=10
            else:    
                call("awk \'{if($NF~/images/ && $1>2) print NR}\' mckq/out'%d'_'%d' > mckq/tp_%d"%(glno,ip,glno),shell=True);
                ## if no sources with no. of lensed images >2, try to find
                ## doubly imaged sources
                if( os.path.getsize("mckq/tp_%d"%(glno))<=0):
                    imgnos=0.4
                    cond=cond+1
                else:
                    rwno=np.loadtxt("mckq/tp_%d"%(glno),unpack=True)
                    cond=10
        if(cond==3):
            print("No lensed images found for lens ",ii,"with id",gid,"..continuing..")
            continue
        else: 
            ## More than one eligible sources found 
            if(rwno.ndim>0):
                 #np.random.seed(myseed)
                 indx=np.random.randint(0,rwno.size)
                 rwn=int(rwno[indx])
            ## One eligible source found 
            else:
                 rwn=int(rwno)
                 
            filename=('mckq/out%d_%d'%(glno,ip))
                
            ## Extract source position and no. of lensed images
            call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d \\n\",$1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s > mckq/tmp_%d"%(filename,glno),shell=1)
            
            ## Extract source position, no. of lensed images, row no. of the
            ## first lensed image and lens galaxy id
            call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d %d '+"%d"%(gid)+' \\n\",$2,NR+1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s >> mckq/srcinfo_%d"%(filename,glno),shell=1)
            xx,yy,jj=np.loadtxt('mckq/tmp_%d'%(glno),unpack=True)
            call("awk \'{if(NR>'%d' && NR<='%d') print $3}\' '%s'  >  mckq/tmp2_%d" % (rwn,rwn+jj,filename,glno),shell=True)
            magn=np.loadtxt('mckq/tmp2_%d'%(glno),unpack=True)
            magnsrt=sorted(np.abs(magn),reverse=1)
            magsum=sum(np.abs(magn))
            
            srcx[ii]=xx
            srcy[ii]=yy
            smag[ii]=magnsrt[1] 
            summ[ii]=magsum 
            imno[ii]=int(jj)
    return srcx,srcy,smag,summ,imno


def getsrcinfo(glno,ip,rwn,gid):                
    #print('Inside getsrcinfo: ',glno,ip,rwn,gid)
    filename=('mckg/out%d_%d'%(glno,ip))
    print("gllens%d.in"%(glno),filename,rwn)
    
    ## Extract source position and no. of lensed images
    call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d \\n\",$1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s > mckg/tmp_%d"%(filename,glno),shell=1);
    ## Extract source position, no. of lensed images, row no. of the
    ## first lensed image and lens galaxy id
    call('awk \'BEGIN{bool=0;}{if(bool){printf(\"%d %d '+"%d"%(gid)+' \\n\",$2,NR+1); bool=0;}'+"if(NR==%d)"%(rwn-1)+'{printf(\"%f %f \",$1,$2); bool=1;}}\''+" %s >> mckg/srcinfo_%d"%(filename,glno),shell=1);
    xx,yy,jj=np.loadtxt('mckg/tmp_%d'%(glno),unpack=True)
    
    call("awk \'{if(NR>'%d' && NR<='%d') print $3}\' mckg/out'%d'_'%d' >  mckg/tmp2_%d" % (rwn,rwn+jj,glno,ip,glno),shell=True);
    magn=np.loadtxt('mckg/tmp2_%d'%(glno),unpack=True)
    magnsrt=sorted(np.abs(magn),reverse=1)
    magsum=sum(np.abs(magn))
    print(xx, yy, int(jj), magnsrt[1], magsum)

    return xx, yy, int(jj), magnsrt[1], magsum 

## Extract source position and some lensed image properties for galaxies
def srcposrng(reinst,ell,ell_pa,sh_str,sh_pa,glno,flag,gid,lenscode):
    limr=reinst.size
    #print('Statring function -- limr',limr,reinst,ell,ell_pa,sh_str,sh_pa,glno,flag,gid,num,'inputs')
    if(flag): #a
        fp=open('mckg/gllens%d.in'%(glno),'w')
 
       #np.savetxt(fp,['set rscale=15'],fmt='%s');
       #np.savetxt(fp,['set ngrid1=40'],fmt='%s');
       #np.savetxt(fp,['set ngrid2=40'],fmt='%s');
       #np.savetxt(fp,['set verbose=0'],fmt='%s');
       #np.savetxt(fp,['set omitcore=0.02'],fmt='%s \n');
        fp.write("set rscale=15 \n set ngrid1=40 \n set ngrid2=40 \n set verbose=0 \n set omitcore=0.02 \n")
 
        for ii in range(limr):
            #np.savetxt(fp,['startup 1 1'],fmt='%s');
            fp.write("startup 1 1\n")
            fp.write('alpha %f 0 0 %f %f %f %f 0 0 1 \n' %(reinst[ii],ell,ell_pa,sh_str,sh_pa))
            fp.write(" 0 0 0 0 0 0 0 0 0 0 \n")
            #np.savetxt(fp,[' 0 0 0 0 0 0 0 0 0 0'],fmt='%s');
            fp.write("mock1 mckg/out%d_%d 500 2 \n \n" %(glno,ii+1))
 
        #np.savetxt(fp,['quit'],fmt='%s');
        fp.write("quit")
        fp.close()

        print("Running mckg/gllens%d.in"%(glno))
        call('%s mckg/gllens%d.in > mckg/gltmpout_%d'%(lenscode,glno,glno),shell=1)
    
    srcx=np.zeros(limr)#np.arange(limr)*0
    srcy=np.zeros(limr)#np.arange(limr)*0
    smag=np.zeros(limr)#np.arange(limr)*0
    summ=np.zeros(limr)#np.arange(limr)*0
    imno=np.zeros(limr)#np.arange(limr)*0
    
    for ii in range(limr):
        ip=ii+1
        #print "###########starting ",ip
        #np.random.seed(myseed)
        imgnos=np.random.uniform(0,1)
        ## Loop to select a source that has 2 or more lensed images such that a
        ## double is chosen 50% of the times
        cond=1
        while(cond<3):
            if(imgnos<0.5):
                call("awk -f scr2img.awk mckg/out'%d'_'%d' > mckg/tp_%d"%(glno,ip,glno),shell=True)
                ## if no doubly images sources (unlikely to happen), try sources with no. of lensed
                ## images >2
                if( os.path.getsize("mckg/tp_%d"%(glno))<=0):
                    imgnos=0.6
                    cond=cond+1
                else:
                    rwno=np.loadtxt("mckg/tp_%d"%(glno),unpack=True)
                    cond=8
                   #print "coming out of while cond,im=2"
            else:    
                call("awk \'{if($NF~/images/ && $1>2) print NR}\' mckg/out'%d'_'%d' > mckg/tp_%d"%(glno,ip,glno),shell=True)
                print(glno, ip)
                ## if no sources with no. of lensed images >2, try to find
                ## doubly imaged sources
                if( os.path.getsize("mckg/tp_%d"%(glno))<=0):
                    imgnos=0.4
                    cond=cond+1
                else:
                    rwno=np.loadtxt("mckg/tp_%d"%(glno),unpack=True)
                    cond=10
                   #print "coming out of while cond,im>2"
        ########## done while loop 
        if(cond==3):
            print("No lensed images found for lens ",ii,"with id",gid,"..continuing..")
            continue
        else:
            if(rwno.ndim==0):    
                rwn=int(rwno)
#               print("only 1 source found",rwn)
            else:## More than one eligible sources found 
                #np.random.seed(myseed)
                indx=np.random.randint(0,rwno.size)
                rwn=int(rwno[indx])
            srcx[ii],srcy[ii], imno[ii], smag[ii],summ[ii]=getsrcinfo(glno,ip,rwn,gid)
            sys.stdout.flush()
    return srcx,srcy,smag,summ,imno


## Calculate size of the source assuming size-luminosity relation
## using Bernardi et al. 2003 Eqns given in Oguri 2006
def srcsize(mapp,zsrc,pixsc,cosmo):
    ## mapp is in g-band
    Dlum=cosmo.luminosity_distance(zsrc).value  #cc.Dlofz(zsrc)/p.hval
    Mabs =mapp-5*log10(Dlum)-25
    Lum_src=10**(-0.4*(Mabs-5.48))
    Da=cosmo.angular_diameter_distance(zsrc).value*1.0e3 #cc.Daofz(zsrc)/p.hval*1.e3 ## Da is in kpc

    Lrat= Lum_src/10**10.2
    Reff= 10**0.52*Lrat**(2./3.) * 1./(1+zsrc)**2
    ## Half-light radius converting from kpc->radian->arcsec->pix
    #print('Inputs for srcsize: ',mapp,zsrc,pixsc)
    #print('Output: ',(Reff/Da)*(180.*3600/pi/pixsc))
    return (Reff/Da)*(180.*3600/pi/pixsc)
