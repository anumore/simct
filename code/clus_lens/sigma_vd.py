#!/usr/bin/env python
import numpy as np
from math import *
from input_g import *

## Calculate lens velocity dispersion using g,r magnitudes and redshift
## Use absolute magnitudes
def getsigma(mg,mr,zz):
    ## Get the SDSS magnitudes from megacam mags:
    ## g_Mega = g_SDSS - 0.153 (g_SDSS - r_SDSS)
    ## r_Mega = r_SDSS - 0.024 (g_SDSS - r_SDSS)
    ## g_Mega-r_Mega=(g_SDSS-r_SDSS)*(1-0.153+0.024)
    ## (g_Mega-r_Mega)/0.871=(g_SDSS-r_SDSS)
    ## r_SDSS = r_Mega + 0.024*(g_Mega-r_Mega)/0.871 

    Dlum=cc.Dlofz(zz);
    Mabsr=mr-5.0*log10(Dlum)-25.0;
     
    ## Get the SDSS magnitude
    mrsdss=Mabsr+0.024*(mg-mr)/0.871;

    ## Transfer to r' = ^{0.1}r Hubble type E, Frei and Gunn
    mrsdss=mrsdss-0.11;

    ## Assume that the luminosity function evolves such that Mr* decline by 1.5
    ## magnitudes from 0.0 to 1.0, this is similar to the evolution in the B band
    ## found by Faber et al. 2007 from DEEP2 and COMBO-17. This is really an adhoc
    ## prescription but nothing better known to me right now.
    
    mrstar=(-20.44)+(zz-0.1)*1.5;
    LbyLstar=10.0**(-0.4*(mrsdss-mrstar));
    
    ## Parker et al. 2007, Table 1 - Bright sample
    return 142.0*LbyLstar**(1./3.);

