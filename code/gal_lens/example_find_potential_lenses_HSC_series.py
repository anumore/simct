import numpy as np
from astropy.table import Table

from astropy.cosmology import FlatLambdaCDM
#from slsim.Util.mag2errors import get_errors_Poisson




# define the cosmology and the universal constants to be used
cosmo = FlatLambdaCDM(H0=72,Om0=0.26)
constants = {'G':4.2994e-9, 'light_speed':299792.458}  # in km^3 Msun^(-1) year^(-2) and km/s respectively

bkggalcatalog    = "../catalogs/hsc_UDEEP_cosmos_bkggal_mod.csv" 

lenscatalog = '../catalogs/hsc_udeep_frggal_mod.csv'
data = Table.read(lenscatalog,format='csv')
zlens,gg,gr,gi,gz,gy= data['photoz_median'],data['gmag'],data['rmag'],data['imag'],data['zmag'],data['ymag']
gerr,rerr,ierr,zerr,yerr= data['gmag_err'],data['rmag_err'],data['imag_err'],data['zmag_err'],data['ymag_err']
gra,gdec =  data['ra'],data['dec']
gid = data['object_id']


gal_ell,gal_pa = data['ellipticity'],data['position_angle']
gal_stellar_mass = data['stellar_mass']

foreground_table = Table()

foreground_table['ra']  = gra
foreground_table['dec']  = gdec
foreground_table['object_id']  = gid
foreground_table['mag_true_g_lsst'] = gg
foreground_table['mag_true_r_lsst'] = gr
foreground_table['mag_true_i_lsst'] = gi
foreground_table['mag_true_z_lsst'] = gz
foreground_table['mag_true_y_lsst'] = gy
foreground_table['gmag_err'] = gerr
foreground_table['rmag_err'] = rerr
foreground_table['imag_err'] = ierr
foreground_table['zmag_err'] = zerr
foreground_table['ymag_err'] = yerr
foreground_table['redshift'] = zlens
foreground_table['ellipticity_true'] = gal_ell
foreground_table['position_angle_true'] = gal_pa
foreground_table['stellar_mass'] = gal_stellar_mass


# define the keywords for the survey, deflector, and the background source to be used
kwargs_source_gal={'type':'galaxy', 'min_z':0.2, 'max_z':4.0, 'min_mag':23.0, 'max_mag':26.0, 'boost_csect':300, 'mag_limit':26.8,
                     'brightness_lim':19.0,    'tol_mag':0.25, 'tol_z_strict':0.07, 
                       'tol_z_relaxed':2.0, 'tol_color_strict':1.5,
                       'ell_min':0.1, 'ell_max': 0.8, 'PA_min': 0., 'PA_max':180., 'R_Einst_min': 0.5, 'R_Einst_max':3.0}
kwargs_source_qso={'type':'qso', 'min_z':0.2, 'max_z':4.0,'min_mag':23.0, 'max_mag':26.0, 'boost_csect':30000, 'mag_limit':26.8,
                     'brightness_lim':19.0, 'tol_mag':1.0, 'tol_z':0.5, 'tol_z_strict':4.0, 
                       'tol_z_relaxed':40.0, 'tol_color_relaxed':2.5,
                       'ell_min':None, 'ell_max': None,  'PA_min': None, 'PA_max':None, 'R_Einst_min': 0.5, 'R_Einst_max':3.0}
kwargs_deflector={'min_z':0.2,'max_z':1.1, 'min_shear':0.001,'max_shear':0.02,'min_PA_shear':0.0,'max_PA_shear':180.0}

psfpath          = "../new_psf"                               ## Path to psf images
lenscode         = "/home/vibhorenegi/soft/lensmodel"                                                  ## Path to keeton's lensmodel executable
kwargs_survey={'name':'HSC','pixscale':0.168, 'zeropoint':27,'bands':np.array(['g','r','i','z','y'])}
image_size = 101

outdir = './outdir'
fitspath = '/home/vibhorenegi/Codes/lens_simulation/SIMCT_HSC_run_Dec2024/test_boost3000_cosmos_images/inpfits'


#fitspath = '/home/vibhorenegi/Codes/lens_simulation/SIMCT_HSC_run_Dec2024/test_small_boost_SHC_3000/inpfits'
#exptime=np.array([600,600,1200,1200,1200])
#exptime=np.array([25200,25200,50400,68040,68040])  #target depths
exptime=np.array([4200,3960,5880,10620,5580])       #PDR3 paper   
'''# define the keywords for the survey, deflector, and the background source to be used
kwargs_source_gal={'type':'galaxy', 'min_z':0.2, 'max_z':5.0, 'min_mag':23.0, 'max_mag':28.0, 'boost_csect':3, 
                     'ell_min':0.1, 'ell_max': 0.8, 'PA_min': 0., 'PA_max':180., 'R_Einst_min': 0.5, 'R_Einst_max':3.0}
kwargs_source_qso={'type':'qso', 'min_z':0.2, 'max_z':5.0,'min_mag':23.0, 'max_mag':28.0, 'boost_csect':30000, 
                     'ell_min':None, 'ell_max': None,  'PA_min': None, 'PA_max':None, 'R_Einst_min': 0.5, 'R_Einst_max':3.0}
kwargs_deflector={'min_z':0.2,'max_z':2.0, 'min_shear':0.001,'max_shear':0.02,'min_PA_shear':0.0,'max_PA_shear':180.0}
#kwargs_survey={'zeropoint_g':28.51,'zeropoint_r':28.36,'zeropoint_i':28.17,'exptime_g':15,'exptime_r':15,'exptime_i':15}
 #kwargs_deflector={'min_z':0.2,'max_z':1.9, 'min_shear':0.001,'max_shear':0.02,'min_PA_shear':0.0,'max_PA_shear':180.0}
'''

deflector_type = "galaxy"
background_source = "galaxy"
#background_source = "quasar"



print('Original size of foreground table: ',len(foreground_table))
'''id_redshift = np.where(foreground_table['redshift']<2)[0]
foreground_table = foreground_table[id_redshift]
print('Size of foreground table: ',len(foreground_table))'''

foreground_table.rename_column('mag_true_g_lsst', 'gmag')
foreground_table.rename_column('mag_true_r_lsst', 'rmag')
foreground_table.rename_column('mag_true_i_lsst', 'imag')
foreground_table.rename_column('mag_true_z_lsst', 'zmag')
foreground_table.rename_column('mag_true_y_lsst', 'ymag')
foreground_table.rename_column('ellipticity_true', 'ellipticity')
#foreground_table.info()

if 'shear' and 'shear PA' in foreground_table.colnames:
    foreground_table.rename_column('shear', 'random shear')
    foreground_table.rename_column('shear PA', 'random shear PA')
    
else:
    # generate random shear to the galaxies if shear is not already provided
    myseed = 2894 #9999
    np.random.seed(myseed)
    random_shear_value=np.random.uniform(low=kwargs_deflector.get('min_shear'),high=kwargs_deflector.get('max_shear'))
    #print('random_shear_value is :',random_shear_value)
    #random_shear = np.full(len(foreground_table),random_shear_value)
    random_shear=np.random.uniform(low=kwargs_deflector.get('min_shear'),high=kwargs_deflector.get('max_shear'),size=len(foreground_table))

    np.random.seed(myseed)
    random_PA_value=np.random.uniform(kwargs_deflector.get('min_PA_shear'),kwargs_deflector.get('max_PA_shear'))
    #random_shear_PA = np.full(len(foreground_table),random_PA_value)
    random_shear_PA = np.random.uniform(kwargs_deflector.get('min_PA_shear'),kwargs_deflector.get('max_PA_shear'),size=len(foreground_table))

    # add the shear and shear PA to the galaxy catalog table
    foreground_table['random shear'] = random_shear
    foreground_table['random shear PA'] = random_shear_PA




foreground_table['err_g'] = foreground_table['gmag_err']
foreground_table['err_r'] = foreground_table['rmag_err']
foreground_table['err_i'] = foreground_table['imag_err']
foreground_table['err_z'] = foreground_table['zmag_err']
foreground_table['err_y'] = foreground_table['ymag_err']

"""
To run the code in serial
--Run it as python example_find_potential_lenses.py

To run the code in parallel
--The code should be run as 'mpirun -np 4 python example_find_potential_lenses.py'
"""

# Run the code in series
run = 'series'
#run = 'parallel'




'''
#np.random.seed(12345)
#ids = np.random.randint(1,len(foreground_table),500)
ids = np.arange(0,len(foreground_table))

## To run in series
if(run=='series'):
    from find_lenses import find_potential_lenses
    find_potential_lenses(foreground_table[ids], background_source='galaxy', 
                        kwargs_source=kwargs_source_gal, kwargs_deflector=kwargs_deflector, 
                        cosmo=cosmo, constants=constants)

    #find_potential_lenses(foreground_table[0:10000], background_source="quasar",
    #                      kwargs_source=kwargs_source_qso,kwargs_deflector=kwargs_deflector,
    #                      cosmo=cosmo,constants=constants)

elif(run=='parallel'):
    from find_lenses import find_potential_lenses_parallel
    from mpi4py import MPI
    ## MPI Initialization
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    find_potential_lenses_parallel(comm,rank,size,foreground_table[ids],background_source="galaxy",kwargs_source=kwargs_source_gal,kwargs_deflector=kwargs_deflector,cosmo=cosmo,constants=constants)

else:
    print('Invalid input')

'''

##############################


'''
sorted_gal_table = Table.read('sorted_galaxies.csv', format='csv')
all_data_table = Table.read('all_sources.csv', format='csv')
ids = np.random.randint(1,len(sorted_gal_table),2000)
sorted_gal_table = sorted_gal_table[ids]#[4170:4180]#[ids]
#print(sorted_gal_table)

#run = 'series'
if(run=='series'):
    from find_lensed_sources import find_lensed_sources_properties
    find_lensed_sources_properties(deflector_type,source_type='galaxy', sorted_gal_table= sorted_gal_table, 
                                   all_data_table=all_data_table, bkgcatalog=bkggalcatalog, kwargs_source=kwargs_source_gal, 
                                   kwargs_deflector=kwargs_deflector, kwargs_survey=kwargs_survey,
                                   lenscode=lenscode,psfpath=psfpath,image_size=image_size,cosmo=cosmo)

                                   
elif(run=='parallel'):
    from find_lensed_sources import find_lensed_sources_properties_parallel
    from mpi4py import MPI
    ## MPI Initialization
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    find_lensed_sources_properties_parallel(comm,rank,size,deflector_type,source_type='galaxy', sorted_gal_table= sorted_gal_table, 
                                   all_data_table=all_data_table, bkgcatalog=bkggalcatalog, kwargs_source=kwargs_source_gal, 
                                   kwargs_deflector=kwargs_deflector, kwargs_survey=kwargs_survey,
                                   lenscode=lenscode,psfpath=psfpath,image_size=image_size,cosmo=cosmo)



'''

'''
if(run=='series'):
    from run_gravlens import rungl
    rungl(lenscode,filename='finalpar.csv')
'''

'''
import os
import shutil

out_dir = './outdir'
if os.path.exists(out_dir):
    shutil.rmtree(out_dir)
else:
    pass
os.makedirs(out_dir)


if(run=='series'):
    import os
    from add_noise import add_Poisson_noise
    os.system('mv LOG*.dat gout/')
    add_Poisson_noise(filename='finalpar.csv',kwargs_source=kwargs_source_gal,kwargs_deflector=kwargs_deflector,
                    kwargs_survey=kwargs_survey,outdir=outdir,exptime=exptime)


'''


'''
if(run=='series'):
    import os
    from get_HSC import get_data
    get_data(filename='finalpar.csv')
    os.system('sh hsc_download.sh')

'''



if(run=='series'):
    from merge_images import addto_realimage
    addto_realimage(filename='finalpar.csv',kwargs_source=kwargs_source_gal,kwargs_deflector=kwargs_deflector,
                    kwargs_survey=kwargs_survey, fitsdir=fitspath, outdir=outdir, exptime=exptime)
    

