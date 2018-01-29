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


# define parameters

diff_limit = 0.06


def func_comp(x, a, b, c):
    x_new = 10**x/c
    return(np.log10(a*x_new**b*np.exp(x_new)))


def get_completeness_correction_numerical(log_Mh, redshift_in):
    '''
    This function returns the completeness correction for
    a given halo mass (all in log units), using the Schechter
    description.
    '''
    redshift_list = np.array([4.0, 6.0, 8.0, 10.0])
    idx_redshift = np.where(redshift_list == redshift_in)[0][0]
    mass_list = completeness_table[0]
    correction_list = completeness_table[idx_redshift+1]
    mass_list = mass_list[~np.isnan(correction_list)]
    correction_list = correction_list[~np.isnan(correction_list)]
    diff = np.interp(log_Mh, mass_list, correction_list, left=correction_list[0], right=correction_list[-1])
    diff[diff < diff_limit] = diff_limit
    return(diff-diff_limit)


def get_completeness_correction_parametrized(log_Mh, redshift_in):
    '''
    This function returns the completeness correction for
    a given halo mass (all in log units), using the Schechter
    description.
    '''
    redshift_list = np.array([4.0, 6.0, 8.0, 10.0])
    ii_redshift = np.abs(redshift_list-redshift_in).argmin()
    if (redshift_list[ii_redshift] == 4.0):
        popt = [1.20813739e+00, 3.04174445e-03, 1.00000000e+13]
    elif (redshift_list[ii_redshift] == 6.0):
        popt = [1.95393928e+00, 8.28309118e-02, 9.99999987e+12]
    elif (redshift_list[ii_redshift] == 8.0):
        popt = [1.47856861e+00, 6.14304840e-02, 2.54652583e+11]
    elif (redshift_list[ii_redshift] == 10.0):
        popt = [1.03167560e+00, 8.87468518e-30, 5.42649366e+10]
    diff = func_comp(log_Mh, popt[0], popt[1], popt[2])
    diff[diff < diff_limit] = diff_limit
    return(diff-diff_limit)



