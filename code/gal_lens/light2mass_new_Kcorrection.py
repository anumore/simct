#!/usr/bin/env python
import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy
from astropy.modeling.models import Linear1D
#from slsim.Util.k_correction import kcorr_sdss
from astropy.cosmology import FlatLambdaCDM
import kcorrect.kcorrect
from color_transformations import HSC_to_DES
from color_transformations import DES_to_SDSS


"""
This module provides function to calculate the central stellar velocity dispersion of the deflector 
(elliptical galaxies) using LSST broadband magnitudes and the redshift. It assumes the evolution of 
galaxy luminosity function, as discussed in Bell et al 2004, Blanton et al 2003, and uses the 
scaling relations for L-sigma relationship from Choi et al 2007 (derived from spectroscopic 
measurements) and Parker et al 2007 (derived from weak lensing measurements). The user has the 
option to decide which scaling relations he/she wants to use, the one derived from spectroscopic
measurements or from the weak lensing measurements.
"""

def kcorr_sdss(mags_sdss,
               redshift,
               abcorrect=False,
               responses=["sdss_u0", "sdss_g0", "sdss_r0", "sdss_i0", "sdss_z0"],
               responses_out=["sdss_u0", "sdss_g0", "sdss_r0", "sdss_i0", "sdss_z0"],
               band_shift=0.0,
               redshift_range=[0, 2],
               cosmo=FlatLambdaCDM(H0=72, Om0=0.26)):
    
    # Etract the magnitudes and errors in separate arrays.
    mags = unumpy.nominal_values(mags_sdss).T
    mag_errs = unumpy.std_devs(mags_sdss).T

    #print('mags',mags.min(),mags.max())
    #print('mag_errs',mag_errs.min(),mag_errs.max())

    responses = ['sdss_g0', 'sdss_r0', 'sdss_i0', 'sdss_z0']
    kc = kcorrect.kcorrect.Kcorrect(responses=responses)

    maggies_ivar = np.zeros(mag_errs.shape, dtype=np.float32)
    maggies = np.zeros(mags.shape, dtype=np.float32)

    mag_low = mags - mag_errs
    mag_high = mags + mag_errs

    #print(mag_low.min(),mag_low.max(),mag_high.min(),mag_high.max())
    for j in range(len(maggies)):
        for k in np.arange(len(responses), dtype=int):
            maggies[j, k] = 10**(-0.4*mags[j, k])
            maggies_ivar[j, k] = 0.5*( 10**(-0.4*mag_low[j, k]) - 10**(-0.4*mag_high[j, k]) )

    # "coeffs" is a [5]-array with coefficients multiplying each template
    coeffs = kc.fit_coeffs(redshift=redshift, maggies=maggies, ivar=maggies_ivar)

    # "k" is a [5]-array with the K-corrections in magnitude units
    k = kc.kcorrect(redshift=redshift, coeffs=coeffs)

    return k




def Lsigma_relation_spectroscopic(mgSDSS, mrSDSS, Dlum, redshift):
    """
    input params:

    mgSDSS: k-corrected g-band magnitude of the deflector
    type:   a 1D array of floats

    mrSDSS: k-corrected r-band magnitude of the deflector
    type:   a 1D array of floats

    Dlum: distance luminosity of the deflector
    type: a 1D array of floats

    redshift: redshift of the deflector
    type: a 1D array of floats

    .. [1] Bell et al., (2004), astro-ph/0303394, doi: 10.1086/420778
    .. [2] Choi, Park and Vogeley, (2007), astro-ph/0611607, doi:10.1086/511060
    """

    # Use the SDSS g-band and r-band magnitudes to get the B-band apparent magnitude of the galaxy using the relation
    # given in equation A2, Appendix, Bell et al 2004 for red galaxies. This is required only for using the relations
    # based on spectroscopic measurements.
    MabsB = mgSDSS + 0.155 + 0.370 * (mgSDSS - mrSDSS)

    # Convert the apparent B-band magnitude to the absolute B-band magnitude using the redshift and cosmology defined
    # Note that the 25 here comes since Dlum is in Mpc
    MabsB = MabsB - 5.0 * np.log10(Dlum/10)# - 25.0
    
    """Now using the data from DEEP2 and COMBO-17 surveys, Bell et 2004 found that the
    B-band luminosity function evolves such that characteristic magnitude MBstar decline
    by 1.5 magnitudes from z=0.0 to z=1.0. We use the same assumption here;

    Hence, MBstar and redshift should follow the relation, i.e., MBstar =
    MBstar0-(redshift)*1.5. where MBstar0 = MBstar(at redshift=0). In our case, MBstar0
    = -19.31 has been estimated from the mean value of the MBstar0, from Table 1, Bell
    et al 2004.
    """
    # define a 1D line model for MBstar evolution with redshift.
    MBstar_func = Linear1D(-1.5, -19.31)

    # Use the above model to calculate MB* at the deflector redshift
    MBstar = MBstar_func(redshift)

    # Calculate L/L* using the magnitude-luminosity relation
    LbyLstar = 10.0 ** (-0.4 * (MabsB - MBstar))
    """Now use the L-sigma relation for the elliptical galaxies i.e., the Faber Jackson
    relation, sigma/sigma_star = (L/Lstar)**(1/alpha) and taking the sigma* and alpha
    value from Choi et al 2007, derived for early type galaxies, calculate the the
    velocity dispersion sigma."""
    sigma_star, alpha = ufloat(161, 5), 2.32  # Choi et al 2007

    # Use sigma_star and alpha values to calculate the stellar velocity dispersion sigma
    sigma = sigma_star * LbyLstar ** (1 / alpha)

    return sigma,MBstar,LbyLstar,MabsB


def Lsigma_relation_weaklensing(mrSDSS, miSDSS, Dlum, redshift):
    """
    input params:

    mrSDSS: k-corrected r-band magnitude of the deflector
    type:   a 1D array of floats

    miSDSS: k-corrected i-band magnitude of the deflector
    type:   a 1D array of floats

    Dlum: distance luminosity of the deflector
    type: a 1D array of floats

    redshift: redshift of the deflector
    type: a 1D array of floats

    .. [1] Blanton et al., (2003), astro-ph/0210215, doi: 10.1086/375776
    .. [2] Parker et al., (2007),  arXiv:0707.1698, doi: 10.1086/521541
    """

    # Convert the apparent r-band magnitudes to the absolute r
    #Mabsr = mrSDSS - 5.0 * np.log10(Dlum/10)# - 25.0

    # Convert the sdss r-mag to r'-mag from Frei & Gunn 2003 (Table 3).
    # r' is a fake filter i.e., r shifted to z=0.1.
    #Mabsr = Mabsr - 0.11
    Mabsr = mrSDSS - 0.11
    """We assume the same assumption here (from Bell et al 2004) for decline of
    characteristic magnitude Mrstar for r'-band,

    Hence, Mrstar and redshift should follow the relation, i.e., Mrstar =
    Mrstar0-(redshift-0.1)*1.5. where Mrstar0 = Mrstar(at redshift=0.1).

    In our case, Mrstar0 = -20.44 has been estimated from Table 2, Blanton et al 2003.
    """

    Mrstar0 = -20.44  # calculated at redhift=0.1
    # Use the above value to calculate Mrstar at the deflector redshift
    Mrstar = Mrstar0 - (redshift - 0.1) * 1.5

    # Calculate L/L* using the magnitude-luminosity relation
    LbyLstar = 10.0 ** (-0.4 * (Mabsr - Mrstar))
    """
    Now use the L-sigma relation and taking the sigma_star and alpha value from
    Parker et al 2007, derived using weak-lensing measurements, calculate the the
    velocity dispersion sigma.
    """
    # sigma_star, alpha = 142+-18, 3      # Parker et al 2007

    sigma_star_nominal = np.ones(len(LbyLstar)) * 142
    sigma_star_stdev = np.ones(len(LbyLstar)) * 18
    alpha = np.ones(len(LbyLstar)) * 3
    sigma_star = unumpy.uarray(sigma_star_nominal, sigma_star_stdev)
    sigma_star[miSDSS > 20.5] = ufloat(137, 11)
    alpha[miSDSS > 20.5] = 3

    # Use sigma_star and alpha values to calculate the stellar velocity dispersion sigma
    sigma = sigma_star * LbyLstar ** (1 / alpha)

    return sigma,Mrstar,LbyLstar,Mrstar


def get_velocity_dispersion(
    deflector_type,
    lsst_mags,
    lsst_errs,
    redshift,
    cosmo=FlatLambdaCDM(H0=70, Om0=0.3),
    bands=["g", "r", "i", "z", "y"],
    scaling_relation="spectroscopic",
):
    """
    input_params:

    deflector_type: type of the foreground/ deflector, e.g., 'elliptical'
    type: string

    lsst_mags: a 2D array of the lsst magnitudes of the deflector with multi-band magnitudes
    along the row, and different deflector along the column.
    type: a 2D array of floats.

    lsst_errs: a 2D array of the lsst magnitude errors of the deflector with multi-band errors
    along the row, and different deflector along the column.
    type: a 2D array of floats.

    Note: Please provide atleast three bands data, including the g, r, and i bands.
    The three bands are required to perform the k-correction in the SDSS bands. If there is some other
    way of doing k-correction directly in the LSST bands, we will need only the two g and r bands data.

    redshift:   a 1D array of the redshifts
    type: a 1D array of floats

    cosmo: cosmology defined
    type: astropy.cosmology

    bands: bands for which you're providing the magnitudes, for now use only 'g', 'r', and 'i'
    type: a list of strings e.g., ['g','r','i']

    returns:
    stellar velocity dispersion [km/s]

    References
    ----------
    .. [1] Bell et al., (2004), astro-ph/0303394, doi: 10.1086/420778
    .. [2] Blanton & Roweis (2007), astro-ph/0606170, doi: 10.1086/510127
            https://kcorrect.readthedocs.io/en/5.1.2/
    .. [3] Choi, Park and Vogeley, (2007), astro-ph/0611607, doi:10.1086/511060
    .. [4] Blanton et al., (2003), astro-ph/0210215, doi: 10.1086/375776
    .. [5] Parker et al., (2007),  arXiv:0707.1698, doi: 10.1086/521541
    """

    if deflector_type != "elliptical":
        raise KeyError("The module currently supports only elliptical galaxies.")

    lsst_bands = ["u", "g", "r", "i", "z", "y"]

    # extract the indices of the available lsst bands
    indices = [lsst_bands.index(band) for band in bands]

    lsst = {}
    for ind in range(len(indices)):
        lsst["{0}".format(lsst_bands[indices[ind]])] = unumpy.uarray(
            lsst_mags[:, ind], lsst_errs[:, ind]
        )

    mgDES,mrDES,miDES,mzDES,myDES = HSC_to_DES(lsst['g'],lsst['r'],lsst['i'],lsst['z'],lsst['y'])
    mgSDSS,mrSDSS,miSDSS,mzSDSS = DES_to_SDSS(mgDES,mrDES,miDES,mzDES,myDES)

    if scaling_relation == "spectroscopic":
        # for k-correction upto redshift z=0 only
        band_shift = 0.0

    elif scaling_relation == "weak-lensing":
        # for k-correction upto redshift z=0.1
        # since the scaling relations used are at z=0.1
        band_shift = 0.1

    else:
        raise KeyError("Invalid input for scaling relations.")

    
    sdss_responses = [ "sdss_g", "sdss_r", "sdss_i", "sdss_%z" ]

    # Find out the K-correction factor using the kcorrect module by Blanton
    k_corrections = kcorr_sdss(
        np.array([mgSDSS, mrSDSS, miSDSS,mzSDSS]),
        redshift,
        responses=sdss_responses,
        responses_out=sdss_responses,
        band_shift=band_shift,
        redshift_range=[0, 2],
    )

    # Apply the K-correction on the SDSS magnitudes
    mgSDSS = mgSDSS - k_corrections[:, 0]
    mrSDSS = mrSDSS - k_corrections[:, 1]
    miSDSS = miSDSS - k_corrections[:, 2]
    mzSDSS = mzSDSS - k_corrections[:, 3]


    ## Note: It will be better if we apply the K-correction directly on the LSST magnitudes,
    ## but no such relation is known to Vibhore right now.

    # calculates the distance luminosity using the redshift and the cosmology
    Dlum = cosmo.luminosity_distance(redshift).to("pc").value
    #Dlum = (cosmo.luminosity_distance(redshift) * cosmo.H(0) / 100).value

    if scaling_relation == "spectroscopic":
        # Use the Lsigma relation based on spectroscopic measurements to calculate the
        # sigma of the deflector
        sigma,Mstar,LbyLstar,Mabs = Lsigma_relation_spectroscopic(mgSDSS, mrSDSS, Dlum, redshift)

    elif scaling_relation == "weak-lensing":
        # Use the Lsigma relation based on weak-lensing measurements to calculate the
        # sigma of the deflector
        sigma,Mstar,LbyLstar,Mabs = Lsigma_relation_weaklensing(mrSDSS, miSDSS, Dlum, redshift)


    #print('sigma,Mstar,LbyLstar', sigma,Mstar,LbyLstar)
    # returns the calculated sigma
    # type: a 1D array of uncertainties.core.Variable
    ##   to extract the nomianl values and the uncertainities in separate arrays,
    ##   use unumpy.nominal_values(sigma) and unumpy.std_devs(sigma)
    return sigma
