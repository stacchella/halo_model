'''
Sandro Tacchella
December 20, 2017
'''

import numpy as np


def read_in_efficency(file_name):
    '''
    Reads in and manipulates dark matter halo accretion history.
    '''
    # read in catalog
    epsilon_info = np.load(file_name)
    epsilon_info_Mh = epsilon_info[0]
    epsilon_info_epsi = epsilon_info[1]
    # define epsilon function
    def draw_epsilon(Mh_in, size_in=1.0):
        '''
        This function returns an efficency from
        the calibrated distribution for a given halo mass.
        '''
        idx_Mh = np.abs(epsilon_info_Mh-Mh_in).argmin()
        epsilon_median = epsilon_info_epsi[idx_Mh]
        return(epsilon_median*np.ones(size_in))
    return(draw_epsilon)


# def read_in_efficency(file_name):
#     '''
#     Reads in and manipulates dark matter halo accretion history.
#     '''
#     # read in catalog

#     epsilon_info = np.load(file_name)
#     epsilon_info_Mh = epsilon_info.T[0][1:]
#     epsilon_info_epsi = epsilon_info[0][1:]
#     epsilon_matrix = epsilon_info[1:, 1:]
#     # define epsilon function
#     def draw_epsilon(Mh_in, size_in=1.0):
#         '''
#         This function returns an efficency from
#         the calibrated distribution for a given halo mass.
#         '''
#         idx_Mh = np.abs(epsilon_info_Mh-Mh_in).argmin()
#         epsilon_distribution = epsilon_matrix[idx_Mh]
#         return(np.random.choice(epsilon_info_epsi, size=size_in, replace=True, p=epsilon_distribution/np.sum(epsilon_distribution))) 
#     return(draw_epsilon)


