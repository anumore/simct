#!/usr/bin/env python
from math import *
from subprocess import call
import multiprocessing
#from mpi4py import MPI
import time
import numpy as np
import sys
import gensrcpos as scp
import genimg as gn
#from input_params import *
import os
from astropy.table import Table

##########################################################################################
## ***********************    INPUTS:   ***********************                         ##
## Foreground galaxy catalog with lens model properties as output by                    ##
## find_lenses.py and background galaxy catalog which has Redshift and Magnitude          ##
##                                                                                      ##
## ***********************    PURPOSE:  ***********************                         ##
## Use the lens properties to choose a single source per lens, calculate                ##
## source position based on the limits on fluxes of the lensed images, assign           ##
## colors and generate input files which will produce lensed galaxy images              ##
##########################################################################################

## Galaxies as background sources

def find_lensed_sources_properties(deflector_type,source_type, sorted_gal_table, 
                                   all_data_table, bkgcatalog, kwargs_source, kwargs_deflector, 
                                   kwargs_survey,lenscode,psfpath,image_size,cosmo,
                                   parallel=False, chunk_offset=0):
    '''
    ## Read the lens galaxy catalog, real background galaxy color catalog and 
    ## bkg galaxy catalog generated for each potential lens 

    params:

    deflector_type: type of the deflector e.g., galaxy-elliptical
    type: string

    source_type: background source for the lensing. e.g., galaxy, or quasar.
    type: string

    bkgcatalog: catalog of the background sources with 
        source_ra (deg), source_dec (deg), source_redshift, xx, xx, xx,
    type: ascii file

    kwargs_source={'type':'galaxy', 'min_z':1.2, 'max_z':4.0, 'min_mag':23.0, 'max_mag':26.0, 'boost_csect':3, 
                'mag_limit':26.8, 'brightness_lim':19.0,   'tol_mag':1.0, 'tol_z':0.5, 'tol_z_strict':4.0, 
                'tol_z_relaxed':20.0, 'tol_color_strict':1.5, 'tol_color_relaxed':1.5, 'ell_min':0.1, 'ell_max': 0.8,  
                'PA_min': 0, 'PA_max':180.0, 'R_Einst_min': 0.5, 'R_Einst_max':3.0}
        all the background source parameters
    type: dictionary
                
    kwargs_deflector={'min_z':0.2,'max_z':1.1,'min_shear':0.001,'max_shear':0.02,'min_PA_shear':0,'max_PA_shear':180.0}
        all the foreground deflector parameters
    type: dictionary

    kwargs_survey={'name':'HSC','pixscale':0.168, 'zeropoint':27,'bands':np.array(['g','r','i','z','y'])},
        all the survey parameters
    type: dictionary

    lenscode: path to the lenscode software
    type: 'string'

    psfpath: path to the psf files
    type: string

    image_size: size of the final images chunk
    type: int

    cosmo: cosmology defined
    type: class 'astropy.cosmology.flat.FlatLambdaCDM'
    
    '''

    # read the sorted_galaxies.csv file for the deflector ra, dec, redshift, g,r,i,z,y mags, ellipticity, 
    # position angle, shear, shear position angle and the unique galaxy id, and fiber id
    #sorted_gal_table = Table.read('sorted_galaxies.csv', format='csv')
    foreground_ra,foreground_dec,foreground_redshift = sorted_gal_table['ra'],sorted_gal_table['dec'],sorted_gal_table['redshift']
    foreground_gmag,foreground_rmag,foreground_imag = sorted_gal_table['gmag'],sorted_gal_table['rmag'],sorted_gal_table['imag']
    foreground_ellipticity,foreground_position_angle = sorted_gal_table['ellipticity'],sorted_gal_table['position_angle_true']
    foreground_shear,foreground_shear_PA = sorted_gal_table['random shear'],sorted_gal_table['random shear PA']
    foreground_galaxy_id = sorted_gal_table['object_id']
    foreground_id = sorted_gal_table['deflector Id']

    #id_deflector,velocity_dispersion_deflector,Einstein_radius_deflector,smii0,zs0 = gal_table['galaxy_id'],gal_table['velocity dispersion'],gal_table['Einstein radius'],gal_table['source_mag'],gal_table['source_redshift']

    # read the 'all_sources.csv' file to read the ids, velocity dispersion, Einstein radius,
    # and the source i-magnitude and redshift for the potential lensing deflectors 
    #all_data_table = Table.read('all_sources.csv',format='csv')
    id_deflector,velocity_dispersion_deflector,Einstein_radius_deflector,smii0,zs0 = all_data_table['Id'], all_data_table['Velocity dispersion'],all_data_table['Einstein radius'],all_data_table['Source mag'],all_data_table['Source z']
    #print('Yes, all_sources.csv is needed')


    # read the background source catalog for the redshift, g,r,i,z,y magnitudes
    ################################################################################  
    background_table = Table.read(bkgcatalog, format='csv')
    source_redshift,source_ymag,source_gmag,source_rmag,source_imag,source_zmag = background_table['photoz_median'],background_table['ymag'],background_table['gmag'],background_table['rmag'],background_table['imag'],background_table['zmag']

    # convert the Einstein radius to pixels
    Einstein_radius_deflector=Einstein_radius_deflector/kwargs_survey.get('pixscale')

    flxrefg=10.0**(-0.4*(kwargs_source.get('mag_limit')-kwargs_survey.get('zeropoint')))
    flxlimbrt=10.0**(-0.4*(kwargs_source.get('brightness_lim')-kwargs_survey.get('zeropoint')))
    
    ## Set the flagg to 0, if you don't want the gravlens input files in mckg to
    ## be recreated
    flagg=1
    

    # fix a seed
    myseed = 29824 + chunk_offset*10
    np.random.seed(myseed)

    ## For all sources properties
    #uniq_deflectors,unique_ids= np.unique(id_deflector, return_index=True)      # id_deflector comes from all_sources.csv
    # uniq_deflectors are unique in all_sources.csv
    # unique_ids are ids in all_sources.csv
    # uniq_deflectors are same as foreground_id
    #print(uniq_deflectors,unique_ids)
    #print(unique_ids.min(),unique_ids.max())
    # unique ids will range from 0 to len(all_sources.csv)

    #print(len(uniq_deflectors))
    #print(np.min(foreground_id-uniq_deflectors),np.max(foreground_id-uniq_deflectors))

    keywords1 = ['einstein_radius', 'velocity_dispersion', 'source_imag', 'source_redshift', 'source_xpos', 'source_ypos', 'source_magnification', 'image_count']
    
    # create a dictionary to store the arrays for the deflector and source properties,
    # make them zero in the beginning 
    deflector_array = {}
    for name in keywords1:
        #deflector_array[name] = np.zeros(uniq_deflectors.size)
        deflector_array[name] = np.zeros(foreground_id.size)  # same as uniq_deflectors ?

    # stored in a different dictionary, here array size will be different
    keywords2 = ['gmag', 'rmag', 'imag', 'zmag', 'ymag', 'ellipticity', 'position_angle', 'half_light_radius']
    source_array={}
    for name in keywords2:
        source_array[name] = np.zeros(foreground_id.size)
        
   
    ## Set the range within which to match magnitudes and redshift of real galaxies 
    ## in order to extract colors 
    
    sorted_indices = []
    ## This loop is run for each lens in the lens catalog
    for local_ii in range(len(foreground_ra)):
        ii = local_ii + chunk_offset
        ######################################################
        ## PART 1- Select one bkg source from multiple sources 
        ######################################################
        id_x = np.where(id_deflector == foreground_id[local_ii])[0]
        kkl,kkh = id_x[0],id_x[-1]       
        print('************************************************************kkl and kkh are : ',kkl,kkh)
 
        #kkl=unique_ids[ii]
        #if(ii==uniq_deflectors.size-1):
        #    kkh=id_deflector.size
        #else:
        #    kkh=unique_ids[ii+1]
        #print('******',kkl,kkh,len(Einstein_radius_deflector[kkl:kkh]),Einstein_radius_deflector[kkl:kkh]*kwargs_survey.get('pixscale'))
        #print('kkl to kkh has all entries from all_sources.csv for the same deflector')

        ## For each source from all_sources.csv for a given deflector, extract source position, 
        ## flux of the 2nd brightest image, total magnification of the lensed source 
        ## and number of lensed images
        srcxt,srcyt,smagt,sumt,imtno=scp.srcposrng(Einstein_radius_deflector[kkl:kkh],foreground_ellipticity[local_ii],foreground_position_angle[local_ii],foreground_shear[local_ii],foreground_shear_PA[local_ii],ii+1,flagg,int(foreground_galaxy_id[local_ii]),lenscode)
        #srcxt,srcyt,smagt,sumt,imtno=scp.srcposrng(Einstein_radius_deflector[kkl:kkh],foreground_ellipticity[local_ii],foreground_position_angle[local_ii],foreground_shear[local_ii],foreground_shear_PA[local_ii],local_ii+1,flagg,int(foreground_galaxy_id[local_ii]),lenscode)

        cnt0=np.arange(kkh-kkl)*0
        jj=0
        kk=kkl
        for kk in range(kkl,kkh):
            hh=kk-kkl
            flx2=10.0**(-0.4*(smii0[kk]-kwargs_survey.get('zeropoint'))) * smagt[hh]
            flxall=10.0**(-0.4*(smii0[kk]-kwargs_survey.get('zeropoint'))) * sumt[hh]

            ## Accept each source, if the 2nd brightest lensed image and sum of flux of all lensed
            ## images is above the set limits
            if(flx2>flxrefg and flxall<flxlimbrt):
                cnt0[jj]=kk
                jj=jj+1
            else:
                print(f"Source skipped for lens ID {foreground_galaxy_id[local_ii]}: flx2={flx2}, flxrefg={flxrefg}, flxall={flxall}, flxlimbrt={flxlimbrt}", flush=True)



        ## Choose one source randomly from the sources which satisfy the flux limits
        if(jj>0):
            #np.random.seed(myseed)
            qq=np.random.randint(0,jj)
        else:
            print("No source with 2nd brightest image above mlim=",kwargs_source.get('mag_limit'),"for lens id:",int(foreground_galaxy_id[local_ii]), flush=True)
            continue

        print('Deflector veolocity dispersion is : ',velocity_dispersion_deflector[cnt0[qq]], flush=True)       # this prints correct

        nq=cnt0[qq]-kkl
        print('Deflector velocity dispersion is : ',velocity_dispersion_deflector[nq],flush=True)              # this prints wrong values
        print('Deflector source redshift is : ',zs0[cnt0[qq]],flush=True)                                      # zs0 comes from all_sources.csv
        print('Deflector einstein radius is : ',Einstein_radius_deflector[cnt0[qq]], flush=True)                      # also prints wrong values, why ?
        print('Lengths ',len(velocity_dispersion_deflector),len(zs0),len(Einstein_radius_deflector), flush=True)
        deflector_array['source_imag'][local_ii]=smii0[cnt0[qq]]
        deflector_array['source_redshift'][local_ii]=zs0[cnt0[qq]]
        deflector_array['einstein_radius'][local_ii]=Einstein_radius_deflector[cnt0[qq]]
        deflector_array['velocity_dispersion'][local_ii]=velocity_dispersion_deflector[cnt0[qq]]
        deflector_array['source_xpos'][local_ii]=srcxt[nq]
        deflector_array['source_ypos'][local_ii]=srcyt[nq]
        deflector_array['source_magnification'][local_ii]=smagt[nq]
        deflector_array['image_count'][local_ii]=imtno[nq]
        
        ########################################################################
        ## PART 2- Extract gal colors from real gal catalog
        ########################################################################
        ll=0
        indx_n=np.arange(source_gmag.size)*0
        for jj in range(source_gmag.size):
            if (abs(deflector_array['source_imag'][local_ii]-source_imag[jj])<=kwargs_source.get('tol_mag') and abs(deflector_array['source_redshift'][local_ii]-source_redshift[jj])<=kwargs_source.get('tol_z_strict') and (source_gmag[jj]-source_imag[jj])<kwargs_source.get('tol_color_strict')):
                indx_n[ll]=jj
                ll=ll+1
        
        indx_n=indx_n[0:ll]
        #np.random.seed(myseed)
        nn=np.random.randint(0,ll)
        ncnt0=indx_n[nn]

        ## Use gal colors only (not the magnitudes) from the catalog
        source_array['imag'][local_ii]=deflector_array['source_imag'][local_ii]
        source_array['gmag'][local_ii]=deflector_array['source_imag'][local_ii]+source_gmag[ncnt0]-source_imag[ncnt0]
        source_array['rmag'][local_ii]=deflector_array['source_imag'][local_ii]+source_rmag[ncnt0]-source_imag[ncnt0]
        source_array['zmag'][local_ii]=deflector_array['source_imag'][local_ii]+source_zmag[ncnt0]-source_imag[ncnt0]
        ## NOTE reading cfhtls  u mag into source_y and we discard source_y later 
        source_array['ymag'][local_ii]=deflector_array['source_imag'][local_ii]+source_ymag[ncnt0]-source_imag[ncnt0]
 
        ## Generate random ellipticities, PA and size for the background source
        ell_min,ell_max = kwargs_source.get('ell_min'),kwargs_source.get('ell_max')
        #np.random.seed(myseed)
        source_array['ellipticity'][local_ii]=np.random.uniform(ell_min,ell_max)
        #np.random.seed(myseed)
        source_array['position_angle'][local_ii]=np.random.uniform(kwargs_source.get('PA_min'),kwargs_source.get('PA_max'))

        #half_light_radius = np.array(half_light_radius,dtype=float)
        source_array['half_light_radius'][local_ii]=scp.srcsize(source_array['gmag'][local_ii],deflector_array['source_redshift'][local_ii],kwargs_survey.get('pixscale'),cosmo=cosmo)

        ## Save catalog with all the lens+source parameters
        sorted_indices.append(local_ii+chunk_offset)

        ## Generate lensmodel input files

        gn.genimg_gg(deflector_array['einstein_radius'][local_ii],foreground_ellipticity[local_ii],
                     foreground_position_angle[local_ii],foreground_shear[local_ii],foreground_shear_PA[local_ii],
                     source_array['gmag'][local_ii],source_array['rmag'][local_ii],source_array['imag'][local_ii],
                     source_array['zmag'][local_ii],source_array['ymag'][local_ii],deflector_array['source_xpos'][local_ii],
                     deflector_array['source_ypos'][local_ii],source_array['ellipticity'][local_ii],
                     source_array['position_angle'][local_ii],source_array['half_light_radius'][local_ii],
                     kwargs_survey.get('zeropoint'),foreground_galaxy_id[local_ii],lenscode,
                     psfpath,image_size)

    if parallel:
        print('************',len(source_array['gmag']),len(sorted_indices))
        sys.stdout.flush()
        return [foreground_galaxy_id,foreground_ra,foreground_dec,foreground_redshift,foreground_gmag,foreground_rmag,
                foreground_imag,foreground_ellipticity,foreground_position_angle,deflector_array['einstein_radius'],
                deflector_array['velocity_dispersion'],foreground_shear,foreground_shear_PA, deflector_array['source_xpos'],
                deflector_array['source_ypos'],source_array['gmag'],source_array['rmag'],source_array['imag'],source_array['zmag'],
                source_array['ymag'],deflector_array['source_redshift'], deflector_array['source_magnification'],
                deflector_array['image_count'],source_array['ellipticity'], source_array['position_angle'],
                source_array['half_light_radius'], sorted_indices]
        

        
    else:

        print(len(source_array['gmag']),len(deflector_array['source_ypos']),len(deflector_array['source_magnification']))
        final_table = Table([foreground_galaxy_id,foreground_ra,foreground_dec,foreground_redshift,foreground_gmag,foreground_rmag,
                             foreground_imag,foreground_ellipticity,foreground_position_angle,deflector_array['einstein_radius'],
                             deflector_array['velocity_dispersion'],foreground_shear,foreground_shear_PA,
                            deflector_array['source_xpos'],deflector_array['source_ypos'],source_array['gmag'],source_array['rmag'],
                            source_array['imag'],source_array['zmag'],source_array['ymag'],deflector_array['source_redshift'],
                            deflector_array['source_magnification'],deflector_array['image_count'],source_array['ellipticity'],
                            source_array['position_angle'],source_array['half_light_radius']],
                            names=('galaxy_id','galaxy_ra','galaxy_dec','galaxy_z','galaxy_gmag',
                            'galaxy_rmag','galaxy_imag','galaxy_ellipticity','galaxy_PA','galaxy_Einstein_radius',
                            'galaxy_velocity_dispersion','galaxy_shear','galaxy_shearPA','source_x','source_y',
                            'source_gmag','source_rmag','source_imag','source_zmag','source_ymag','source_redshift',
                            'source_magnification', 'image_count','ellipticity','position_angle','half_light_radius'))
        
        final_table = final_table[sorted_indices]
        final_table.write('finalpar.csv',format='csv',overwrite=True)



def find_lensed_sources_properties_parallel(comm, rank, size, deflector_type,source_type, sorted_gal_table, 
                                   all_data_table, bkgcatalog, kwargs_source, kwargs_deflector, 
                                   kwargs_survey,lenscode,psfpath,image_size,cosmo,
                                   parallel=True):

    # Scatter the data to all processes
    sorted_gal_table_chunks = np.array_split(sorted_gal_table, size)
    local_chunk = comm.scatter(sorted_gal_table_chunks, root=0)

    # Calculate the offset for this chunk
    chunk_offset = sum(len(chunk) for chunk in sorted_gal_table_chunks[:rank])

    print(f"Process {rank} received chunk with {len(local_chunk)} elements")

    # Each process works on its chunk
    local_results = find_lensed_sources_properties(deflector_type,source_type, local_chunk, all_data_table, bkgcatalog, kwargs_source, kwargs_deflector, 
                                   kwargs_survey,lenscode,psfpath,image_size,cosmo,parallel=True,chunk_offset=chunk_offset)

    # Gather results from all processes
    all_results = comm.gather(local_results, root=0)

    if rank == 0:
        # Combine the results from all processes
        foreground_galaxy_id,foreground_ra,foreground_dec,foreground_redshift = [],[],[],[]
        foreground_gmag,foreground_rmag,foreground_imag = [],[],[]
        foreground_ellipticity,foreground_position_angle = [], []
        deflector_einstein_radius,deflector_velocity_dispersion = [],[]
        foreground_shear,foreground_shear_PA = [],[] 
        deflector_source_xpos, deflector_source_ypos= [],[]
        source_gmag,source_rmag,source_imag,source_zmag,source_ymag= [],[],[],[],[]
        deflector_source_redshift, deflector_source_magnification = [],[]
        deflector_image_count= []
        source_ellipticity, source_position_angle, source_half_light_radius = [],[],[]

        sorted_indices = []

        print(len(all_results))
        for result in all_results:
            foreground_galaxy_id.extend(result[0])
            foreground_ra.extend(result[1])
            foreground_dec.extend(result[2])
            foreground_redshift.extend(result[3])
            foreground_gmag.extend(result[4])
            foreground_rmag.extend(result[5])
            foreground_imag.extend(result[6])
            foreground_ellipticity.extend(result[7])
            foreground_position_angle.extend(result[8])
            deflector_einstein_radius.extend(result[9]*0.168)
            deflector_velocity_dispersion.extend(result[10])
            foreground_shear.extend(result[11])
            foreground_shear_PA.extend(result[12])
            deflector_source_xpos.extend(result[13])
            deflector_source_ypos.extend(result[14])
            source_gmag.extend(result[15])
            source_rmag.extend(result[16])
            source_imag.extend(result[17])
            source_zmag.extend(result[18])
            source_ymag.extend(result[19])
            deflector_source_redshift.extend(result[20])
            deflector_source_magnification.extend(result[21])
            deflector_image_count.extend(result[22])
            source_ellipticity.extend(result[23])
            source_position_angle.extend(result[24])
            source_half_light_radius.extend(result[25])
            sorted_indices.extend(result[26])

        print(len(foreground_galaxy_id),len(foreground_ra),len(foreground_dec),len(foreground_redshift))
        print(len(foreground_gmag),len(foreground_rmag),len(foreground_imag))
        print(len(foreground_ellipticity),len(foreground_position_angle),len(deflector_einstein_radius))
        print(len(deflector_velocity_dispersion),len(foreground_shear),len(foreground_shear_PA),len(deflector_source_xpos),len(deflector_source_ypos))
        print(len(deflector_source_redshift),len(deflector_source_magnification),len(deflector_image_count))
        print(len(source_ellipticity),len(source_position_angle),len(source_half_light_radius),len(sorted_indices))
        final_table = Table([foreground_galaxy_id,foreground_ra,foreground_dec,foreground_redshift,foreground_gmag,foreground_rmag,
                             foreground_imag,foreground_ellipticity,foreground_position_angle,deflector_einstein_radius,
                             deflector_velocity_dispersion,foreground_shear,foreground_shear_PA,
                            deflector_source_xpos,deflector_source_ypos,source_gmag,source_rmag,
                            source_imag,source_zmag,source_ymag,deflector_source_redshift,
                            deflector_source_magnification,deflector_image_count,source_ellipticity,
                            source_position_angle,source_half_light_radius],
                            names=('galaxy_id','galaxy_ra','galaxy_dec','galaxy_z','galaxy_gmag',
                            'galaxy_rmag','galaxy_imag','galaxy_ellipticity','galaxy_PA','galaxy_Einstein_radius',
                            'galaxy_velocity_dispersion','galaxy_shear','galaxy_shearPA','source_x','source_y',
                            'source_gmag','source_rmag','source_imag','source_zmag','source_ymag','source_redshift',
                            'source_magnification','image_count','ellipticity','position_angle','half_light_radius',
                            ))
        
        final_table = final_table[sorted_indices]
        final_table.write('finalpar.csv',format='csv',overwrite=True)

    return None
