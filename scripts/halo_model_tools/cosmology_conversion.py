# calculate HMF for cosmology conversion
from astropy.cosmology import WMAP7 as cosmo
from astropy.cosmology import Planck15 as cosmoP
from hmf import MassFunction

# construct analytical halo mass function
# 1 = Press-Schechter, 2 = Jenkins, 3 = Sheth-Tormen, 4 = Warren, 5 = Tinker
# output:
# masses_analytic: Masses used for the generation of the analytic mass function. Units: M_solar
# n_cumulative_analytic: Number density of halos with mass greater then the corresponding mass in masses_analytic. Units: comoving Mpc^-3
# dndM_dM_analytic: Differential number density of halos, (dn/dM)*dM.

# define function for construction of halo mass function


def compute_HMF(redshift, hmf_model='SMT'):
    '''
    hmf_model = PS, Jenkins, SMT, Warren, Tinker08
    '''
    HMF_WMAP = MassFunction(cosmo_model=cosmo, z=redshift, Mmin=1, Mmax=13, hmf_model=hmf_model)
    HMF_Planck = MassFunction(cosmo_model=cosmoP, z=redshift, Mmin=1, Mmax=13, hmf_model=hmf_model)
    #hmf.update(hmf_model='Sheth-Tormen') #update baryon density and redshift
    masses = HMF_WMAP.m/cosmo.h
    MF_WMAP = HMF_WMAP.dndlog10m*cosmo.h**3
    MF_Planck = HMF_Planck.dndlog10m*cosmo.h**3
    return(masses, MF_WMAP, MF_Planck)

