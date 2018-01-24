'''
Sandro Tacchella
December 20, 2017
'''

import numpy as np


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


def construct_SFH(mass_growth_list, t_snapshots, look_back=500.0, dt=0.1, SFH_type=None, epsilon_fct=None, descaling_eff_merg=0.2):
    '''
    This function constructs SFH (time since Big Bang in Myr and SFR) for a given mass growth list
    and snapshot times.
    ---
    mass_list   : mass growth history
    t_snapshots : age of the universe of snapshots
    lookback    : how many Myr back in time to construct SFH
    dt          : time spacing
    SFH_type    : type of star formation history in the last, most recent time bins
                  options: constant, random
    epsilon_fct : efficency function that depends on halo mass
    '''
    # snapshots to consider
    idx_z = np.where(t_snapshots >= t_snapshots[0]-look_back)[0]
    # construct age and mass list
    t_age = np.round(t_snapshots[idx_z]-t_snapshots[idx_z[-1]], 1)[::-1]
    M_list = mass_growth_list[idx_z][::-1]
    # construct SFH:
    time_start = np.round(np.min(t_snapshots[idx_z]), 1)
    t_age_list = np.array([-1.0*time_start, -0.5*dt])
    SFR_list = np.array([0.0, 0.0])
    for ii in range(len(idx_z)-1):
        time_bins = np.arange(0, t_age[ii+1]-t_age[ii]+dt, dt)
        time_center = time_bins[:-1]+0.5*np.diff(time_bins)
        t_age_list = np.append(t_age_list, t_age[ii]+time_center)
        delta_M = np.max([0.0, M_list[ii+1]-M_list[ii]])
        # set SFR distribution in recent time
        if (ii == len(idx_z)-2) and (SFH_type is not None):  # last bin should be according to type
            if (SFH_type == 'constant'):
                mu_SFR, sig_SFR = 0.5*(t_age[ii+1]-t_age[ii]), 1000
            elif (SFH_type == 'random'):
                mu_SFR = np.max(time_bins)*np.random.random(1)
                sig_SFR = 20.0*np.random.random(1)
        else:  # all other bins constant SFR
            mu_SFR, sig_SFR = 0.5*(t_age[ii+1]-t_age[ii]), 1000
        if epsilon_fct is None:
            SFR = model_SFH(time_center, delta_M, mu_SFR, sig_SFR)
        else:
            efficency = epsilon_fct(np.log10(0.5*(M_list[ii+1]+M_list[ii])), size_in=1.0)
            # add here term to lower efficency in merging halos, when not doing calibration
            if (M_list[ii+1] > 2.0*M_list[ii]) & (efficency < 0.0):  # halo is merging
                efficency = efficency + np.log10(descaling_eff_merg)      # because in log units
            SFR = 10**efficency*model_SFH(time_center, delta_M, mu_SFR, sig_SFR)
        SFR_list = np.append(SFR_list, SFR)
    # time_list: array of time since Big Bang in Myr
    time_list = time_start + t_age_list
    SFR_list[~np.isfinite(SFR_list)] = 0.0
    return(time_list, SFR_list)

