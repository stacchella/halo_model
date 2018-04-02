'''
Sandro Tacchella
December 20, 2017
April 02, 2018 (update to include contamination flag)
'''

import numpy as np
import h5py
import glob
import os


path_main = os.environ['WDIR_halo_model']
path_DM_cat = path_main + 'catalogs/DM/'


def get_DM_file(redshift):
    '''
    Get DM filename from redshift.
    '''
    DM_cat_name_list = glob.glob(path_DM_cat + '*.hdf5')
    redshift_list = []
    for ii_name in DM_cat_name_list:
        redshift_list = np.append(redshift_list, float(ii_name.split('_')[-2][1:]))
    return(DM_cat_name_list[np.abs(redshift_list-redshift).argmin()])


def read_in_halo_cat(redshift, cosmo_in):
    '''
    Reads in and manipulates dark matter halo accretion history.
    '''
    # read in catalog
    file_name = get_DM_file(redshift)
    file_read = h5py.File(file_name, 'r')   # 'r' means that hdf5 file is open in read-only mode
    z_table_zX = file_read['Subhaloes/redshift'][:]
    M_table_zX = file_read['Subhaloes/nodeMass'][:]
    Mt_table_zX = file_read['Subhaloes/mbpMass'][:]
    is_contam = file_read['Subhaloes/isContam'][:]
    # convert to solar mass units
    M_table_zX = M_table_zX/cosmo_in.h
    Mt_table_zX = Mt_table_zX/cosmo_in.h
    # reshape Mt table
    Mt_table_zX = Mt_table_zX.reshape(len(Mt_table_zX)/(len(z_table_zX)+1), len(z_table_zX)+1)
    # add initial mass to array (and remove 0 as the end)
    Mt_table_zX = np.vstack([M_table_zX, Mt_table_zX[:-1].T]).T
    # ensure same length as redshift
    Mt_table_zX = Mt_table_zX[:, :len(z_table_zX)-Mt_table_zX.shape[1]]
    return(z_table_zX, M_table_zX, Mt_table_zX, is_contam)

