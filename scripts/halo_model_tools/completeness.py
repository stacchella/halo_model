'''
Sandro Tacchella
January 27, 2018
'''

import numpy as np


def func_comp(x, a, b, c):
    x_new = 10**x/c
    return(np.log10(a*x_new**b*np.exp(x_new)))


def get_completeness_correction(log_Mh, redshift_in):
    '''
    This function returns the completeness correction for
    a given halo mass (all in log units).
    '''
    diff_limit = 0.06
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
    return(diff)



