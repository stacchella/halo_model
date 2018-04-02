'''
Sandro Tacchella
January 27, 2018
'''

# import modules

import numpy as np
import os

# load completeness correction table

path_main = os.environ['WDIR_halo_model']
path_DM_cat = path_main + 'catalogs/DM/'
completeness_table = np.load(path_DM_cat + 'completeness.npy')
completeness_param_table = np.load(path_DM_cat + 'completeness_param.npy')
redshift_list = np.array([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0])


# define parameters

#diff_limit = 0.06


def func_comp(x, a, b, c):
    x_new = 10**x/c
    return(np.log10(a*x_new**b*np.exp(x_new)))


def get_completeness_correction_numerical(log_Mh, redshift_in):
    '''
    This function returns the completeness correction for
    a given halo mass (all in log units), using the Schechter
    description.
    '''
    idx_redshift = np.abs(redshift_list-redshift_in).argmin()
    mass_list = completeness_table[0]
    correction_list = completeness_table[idx_redshift+1]
    mass_list = mass_list[~np.isnan(correction_list)]
    correction_list = correction_list[~np.isnan(correction_list)]
    diff = np.interp(log_Mh, mass_list, correction_list, left=correction_list[0], right=correction_list[-1])
    #diff[diff < diff_limit] = diff_limit
    return(diff)


def get_completeness_correction_parametrized(log_Mh, redshift_in):
    '''
    This function returns the completeness correction for
    a given halo mass (all in log units), using the Schechter
    description.
    '''
    idx_redshift = np.abs(redshift_list-redshift_in).argmin()
    popt = completeness_param_table[idx_redshift]
    diff = func_comp(log_Mh, popt[0], popt[1], popt[2])
    #diff[diff < diff_limit] = diff_limit
    return(diff)



