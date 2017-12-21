'''
Sandro Tacchella
December 21, 2017
'''

import numpy as np


def get_UVandHa_lum(sp_in, age_in):
    '''
    Function computes UV luminosities from model sp.
    - sp:     stellar population model
    - age_in: age in Gyr
    '''
    # get UV magnitudes
    UV_mags = sp_in.get_mags(tage=age_in, bands=['i1500', 'i2300', 'i2800'])
    # convert AB mags to luminosity
    flux_1500 = np.power(10, -0.4*(UV_mags[0]+48.6))
    L_1500 = flux_1500*4*np.pi*(3.086e+19)**2
    flux_2300 = np.power(10, -0.4*(UV_mags[1]+48.6))
    L_2300 = flux_2300*4*np.pi*(3.086e+19)**2
    flux_2800 = np.power(10, -0.4*(UV_mags[2]+48.6))
    L_2800 = flux_2800*4*np.pi*(3.086e+19)**2
    # get Halpha luminosity
    idx_Halpha = (np.abs(sp_in.emline_wavelengths-6564.61)).argmin()
    idx_Hbeta = (np.abs(sp_in.emline_wavelengths-4862.68)).argmin()
    L_Ha = 3.839*10**33*sp_in.emline_luminosity[idx_Halpha]
    L_Hb = 3.839*10**33*sp_in.emline_luminosity[idx_Hbeta]
    return(L_1500, L_2300, L_2800, L_Ha, L_Hb)


def get_luminosities_for_SFH(sp_in, SFH_in, age_list):
    '''
    Function computes Ha and UV luminosities for certain SFH.
    - SFH_in    : 2xN array with age and SFR
    - age_list  : list of ages for which UV and Ha luminosities should be computed
                  in Gyr
    '''
    # set SFH
    sp_in.set_tabular_sfh(SFH_in[0], SFH_in[1], Z=None)
    lum_1500, lum_2300, lum_2800, lum_Ha, lum_Hb = [], [], [], [], []
    for ii_age in age_list:
        #print 'progress [in %]: ', np.round(100.0*np.sum(age_list<=ii_age)/len(age_list))
        lum = get_UVandHa_lum(sp_in, 0.999*ii_age)
        lum_1500 = np.append(lum_1500, lum[0])
        lum_2300 = np.append(lum_2300, lum[1])
        lum_2800 = np.append(lum_2800, lum[2])
        lum_Ha = np.append(lum_Ha, lum[3])
        lum_Hb = np.append(lum_Hb, lum[4])
    return(lum_Ha, lum_Hb, lum_1500, lum_2300, lum_2800)

