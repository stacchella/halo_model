'''
Sandro Tacchella
December 20, 2017
'''

import numpy as np
from astropy.cosmology import WMAP7 as cosmo


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


def construct_SFH(mass_growth_list, t_snapshots, SFH_type=None, epsilon_fct=None, dt_high_res=0.1, dt_low_res=20.0, time_delay=0.1, specific_growth_threshold=0.5):
    '''
    This function constructs SFH (time since Big Bang in Myr and SFR) for a given mass growth list
    and snapshot times.
    ---
    mass_list   : mass growth history
    t_snapshots : age of the universe of snapshots
    SFH_type    : type of star formation history in the last, most recent time bins
                  options: constant, random
    epsilon_fct : efficency function that depends on halo mass
    dt_low_res  : time spacing at early times
    dt_high_res : time spacing at late times
    time_delay  : time delay of DM accretion to baryons in galaxy (fraction of t_H)
    specific_growth_threshold  : threshold of maximum accretion
    '''
    # make sure that mass and time increases, make time bins
    time_bins = np.append(0.0, t_snapshots[:-1][::-1])
    time_center = time_bins[:-1]+0.5*np.diff(time_bins)
    M_growth = mass_growth_list[::-1]
    dM = np.diff(M_growth)
    # specific growth
    dM_M = dM/M_growth[:-1]
    # ensure only positive growth and limit specifc growth
    dM_final = dM.copy()
    dM_final[dM_final < 0.0] = 0.0
    dM_final[dM_M >= specific_growth_threshold] = specific_growth_threshold*M_growth[:-1][dM_M >= specific_growth_threshold]
    # iterate over time
    time_high_resolution = []
    SFR_list = []  # in Msun/yr
    for ii_bin in range(len(time_center)):
        # low resolution regime
        if (time_bins[-1]-time_center[ii_bin] > 200.0):
            time_now = np.arange(time_bins[ii_bin], time_bins[ii_bin+1], dt_low_res)
        # high resolution regime
        else:
            time_now = np.arange(time_bins[ii_bin], time_bins[ii_bin+1], dt_high_res)
        time_high_resolution = np.append(time_high_resolution, time_now)
        SFR_list = np.append(SFR_list, 10**epsilon_fct(np.log10(M_growth[ii_bin]))*cosmo.Ob0/cosmo.Om0*dM_final[ii_bin]/(10**6*(time_bins[ii_bin+1]-time_bins[ii_bin]))*np.ones(len(time_now)))
    # add time delay
    time_vector_shift = time_delay*time_high_resolution
    SFR_list_shifted = np.interp(time_high_resolution, time_high_resolution+time_vector_shift, SFR_list, left=0.0, right=0.0)
    # add special featured SFR
    if (SFH_type == 'random'):
        # re-distribute accretion rate in past 15 Myr
        idx_recent_past = (time_high_resolution > time_high_resolution[-1]-15)
        mass_formed = np.trapz(SFR_list_shifted[idx_recent_past], 10**6*time_high_resolution[idx_recent_past])
        mu_SFR = np.random.choice(time_high_resolution[idx_recent_past])
        sig_SFR = 20.0*np.random.random(1)
        SFR = model_SFH(time_high_resolution, mass_formed, mu_SFR, sig_SFR)
        SFR_final = SFR_list_shifted.copy()
        SFR_final[idx_recent_past] = 0.0
        SFR_final += SFR
    elif (SFH_type == 'constant'):
        SFR_final = SFR_list_shifted.copy()
    SFR_final[~np.isfinite(SFR_final)] = 0.0
    return(time_high_resolution, SFR_final)


