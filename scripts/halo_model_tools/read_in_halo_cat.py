'''
Sandro Tacchella
December 20, 2017
'''

import numpy as np
import h5py


def read_in_halo_cat(file_name, cosmo_in):
    '''
    Reads in and manipulates dark matter halo accretion history.
    '''
    # read in catalog
    file_read = h5py.File(file_name, 'r')   # 'r' means that hdf5 file is open in read-only mode
    z_table_zX = file_read['Subhaloes/redshift'][:]
    M_table_zX = file_read['Subhaloes/nodeMass'][:]
    Mt_table_zX = file_read['Subhaloes/mbpMass'][:]
    # convert to solar mass units
    M_table_zX = M_table_zX/cosmo_in.h
    Mt_table_zX = Mt_table_zX/cosmo_in.h
    # reshape Mt table
    Mt_table_zX = Mt_table_zX.reshape(len(Mt_table_zX)/(len(z_table_zX)+1), len(z_table_zX)+1)
    # add initial mass to array (and remove 0 as the end)
    Mt_table_zX = np.vstack([M_table_zX, Mt_table_zX[:-1].T]).T
    # ensure same length as redshift
    Mt_table_zX = Mt_table_zX[:, :len(z_table_zX)-Mt_table_zX.shape[1]]
    return(z_table_zX, M_table_zX, Mt_table_zX)

