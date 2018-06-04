'''
Sandro Tacchella
December 20, 2017
'''

import numpy as np
from astropy.cosmology import WMAP7 as cosmo
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.cosmology import z_at_value
import astropy.units as u


def gaussian(x, mu, sig):
    return 1.0/(np.sqrt(2.0*np.pi)*sig)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def model_SFH(time, mass_tot, mu, sig):
    '''
    time:      in Myr
    mass_tot:  mass that must be produced
    mu:        epoch when most stars form (in Myr)
    sig:       spread in Myr
    '''
    # get normalization
    norm = mass_tot/np.trapz(gaussian(time, mu, sig), time)
    SFR_list = norm*gaussian(time, mu, sig)*10**-6  # in Msun per yr
    return(SFR_list)


def Mgas_Mstar_relation(redshift, Mstar):
    '''
    Based on Tacconi+18
    '''
    A = np.random.normal(loc=-1.25, scale=0.03)
    B = np.random.normal(loc=2.6, scale=0.25)
    D = np.random.normal(loc=-0.36, scale=0.03)
    log_MgasMstar = A + B*np.log10(1.0+redshift) + D*(np.log10(Mstar)-0.3-10.7)
    return(np.min([3.0, log_MgasMstar]))


def get_dZ(SFR, dmgas_dt, Z, Z0, R, y, lam):
    '''
    Basic bathtub equation a la Lilly (non-equilirium).
    '''
    return((y*(1.0-R)-Z*(1.0-R+lam))*SFR+Z0*dmgas_dt)


def construct_SFH(mass_growth_list, t_snapshots, SFH_type=None, epsilon_fct=None, dt_high_res=0.1, dt_low_res=20.0, time_delay=0.1, time_smoothing=0.02, specific_growth_threshold=1.0, Z0=0.0143*10**-3, R=0.1, y=0.023, lam10=0.4):
    '''
    This function constructs SFH (time since Big Bang in Myr and SFR) for a given mass growth list
    and snapshot times.
    ---
    mass_list   : mass growth history
    t_snapshots : age of the universe of snapshots
    SFH_type    : type of star formation history in the last, most recent time bins
                  options: constant, random
                  => currently not implemented
    epsilon_fct : efficency function that depends on halo mass
    dt_low_res  : time spacing at early times
                  => currently not implemented
    dt_high_res : time spacing at late times
                  => currently not implemented
    time_delay  : time delay of DM accretion to baryons in galaxy (fraction of t_H)
    specific_growth_threshold  : threshold of maximum accretion
    '''
    # make sure that mass and time increases, make time bins
    time_bins = np.append(0.0, t_snapshots[:-1][::-1])
    time_center = time_bins[:-1]+0.5*np.diff(time_bins)
    z_center = np.array([z_at_value(cosmo.age, age*u.Myr) for age in time_center])
    M_growth = mass_growth_list[::-1]
    dM = np.diff(M_growth)
    # specific growth
    dM_M = dM/M_growth[:-1]
    # ensure only positive growth and limit specifc growth
    dM_final = dM.copy()
    dM_final[dM_final < 0.0] = 0.0
    dM_final[dM_M >= specific_growth_threshold] = specific_growth_threshold*M_growth[:-1][dM_M >= specific_growth_threshold]
    # smooth dM_final
    dM_final = convolve(dM_final, Gaussian1DKernel(stddev=time_smoothing*time_center[-1]/np.diff(time_center)[-1]))
    # iterate over time
    epsilon_list = []
    SFR_list = []  # in Msun/yr
    for ii_bin in range(len(time_center)):
        dmgas_dt = cosmo.Ob0/cosmo.Om0*dM_final[ii_bin]/(10**6*(time_bins[ii_bin+1]-time_bins[ii_bin]))
        epsilon_list = np.append(epsilon_list, 10**epsilon_fct(np.log10(M_growth[ii_bin])))
        SFR_list = np.append(SFR_list, epsilon_list[-1]*dmgas_dt)
    # add time delay
    time_vector_shift = time_delay*time_center
    SFR_list_shifted = np.interp(time_center, time_center+time_vector_shift, SFR_list, left=0.0, right=0.0)
    # add special featured SFR
    SFR_final = SFR_list_shifted.copy()
    SFR_final[~np.isfinite(SFR_final)] = 0.0
    # compute Z evolution
    mZ_list = np.array([10**0])
    SFR = np.interp(time_bins, time_center, SFR_final, left=SFR_final[0], right=SFR_final[-1])
    epsilon_list[~np.isfinite(epsilon_list) | (epsilon_list <= 0.0)] = 10**-4
    epsilon_list = np.append(10**-4, epsilon_list)
    dmgas_dt = SFR/epsilon_list
    Z_list = np.array([Z0])
    for idx in range(len(time_bins)-1):
        mstar = 10**3+np.trapz(SFR_final[:idx+1], 10**6*time_bins[:idx+1])
        mgas = 10**Mgas_Mstar_relation(z_center[idx], mstar)*mstar
        lam = lam10*(mstar/10**10)**(-0.33)
        dt = 10**6*(time_bins[idx+1]-time_bins[idx])
        dZ = dt*get_dZ(SFR[idx], dmgas_dt[idx], Z_list[-1], Z0, R, y, lam)
        if np.isnan(dZ):
            mZ_list = np.append(mZ_list, mZ_list[-1])
        else:
            mZ_list = np.append(mZ_list, mZ_list[-1]+dZ)
    Z_list = np.append(Z_list, mZ_list[-1]/mgas)
    return(time_center, SFR_final, mZ_list[1:], Z_list[1:])



