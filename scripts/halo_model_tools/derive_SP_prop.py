'''
Sandro Tacchella
December 21, 2017
'''

import numpy as np

def get_luminosities_for_SFH(sp_in, SFH_in, age_in, idx_EL, wave_interpolate):
    '''
    Function computes Ha and UV luminosities for certain SFH.
    - SFH_in    : 2xN array with age and SFR
    - age_in    : age where compute spectrum, in Gyr
    '''
    # set SFH
    sp_in.set_tabular_sfh(SFH_in[0], SFH_in[1], Z=None)
    # get mass live
    stellar_mass = sp_in.stellar_mass
    # get magnitudes from spectraum WITH emission lines
    sp_in.params['nebemlineinspec'] = True
    mag_list = sp_in.get_mags(tage=age_in, bands=['i1500', 'i2300', 'i2800', 'v', 'u', '2mass_j'])
    L_filters_list = np.power(10, -0.4*(mag_list+48.6))*4*np.pi*(3.086e+19)**2
    # get emission lines
    L_EL_list = 3.839*10**33*sp_in.emline_luminosity[idx_EL]
    # downgrade spectrum without emission lines
    sp_in.params['nebemlineinspec'] = False
    wave, spec = sp_in.get_spectrum(tage=age_in)
    spec_interpolate = np.interp(wave_interpolate, wave, spec)
    return(stellar_mass, L_filters_list, L_EL_list, spec_interpolate)

