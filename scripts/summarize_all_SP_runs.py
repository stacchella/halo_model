'''
Sandro Tacchella
January 3, 2018   : iniate
run this script after run_derive_SP.py in order to
combine all files to one file.

'''

# import modules

import numpy as np
import argparse
import os
import h5py


# define paths

path_main = os.environ['WDIR_halo_model']
path_figures = path_main + 'Figures/Tests/'
path_DM_cat = path_main + 'catalogs/DM/'
path_SFH_cat = path_main + 'catalogs/SFH/'
path_SP_cat = path_main + 'catalogs/SP/'


# set parameters

SFH_type_option = 'constant'  # 'constant' or 'random'

# z=4
# filename_SFH_file = 'SFH_z4_' + SFH_type_option + '.hdf5'
# filename_SP_file = 'SFH_z4_' + SFH_type_option + '_with_L.hdf5'
# filename_SFH_file = 'SFH_z4_' + SFH_type_option + '_calibration.hdf5'
# filename_SP_file = 'SFH_z4_' + SFH_type_option + '_calibration_with_L.hdf5'
# z=6
filename_SFH_file = 'SFH_z6_' + SFH_type_option + '.hdf5'
filename_SP_file = 'SFH_z6_' + SFH_type_option + '_with_L.hdf5'
# z=8
# filename_SFH_file = 'SFH_z8_' + SFH_type_option + '.hdf5'
# filename_SP_file = 'SFH_z8_' + SFH_type_option + '_with_L.hdf5'
# z=10
# filename_SFH_file = 'SFH_z10_' + SFH_type_option + '.hdf5'
# filename_SP_file = 'SFH_z10_' + SFH_type_option + '_with_L.hdf5'


# get number of bins

parser = argparse.ArgumentParser()
parser.add_argument("--number_of_bins", type=int, help="number of cores")
args = parser.parse_args()

number_of_bins = args.number_of_bins


# read in SFH

SFH_file = h5py.File(path_SFH_cat + filename_SFH_file, 'r')
SFH_time = SFH_file['SFH/SFH_time'][:]
SFH_SFR = SFH_file['SFH/SFH_SFR'][:]


# set up new hdf5 file for saving luminosities optained from SFH

try:
    os.remove(path_SP_cat + filename_SP_file)
except OSError:
    pass


lum_file = h5py.File(path_SP_cat + filename_SP_file, 'w')
# add SFH
grp_SFH = lum_file.create_group("SFH")
grp_SFH.create_dataset('SFH_time', data=SFH_file['SFH/SFH_time'][:])
grp_SFH.create_dataset('SFH_SFR', data=SFH_file['SFH/SFH_SFR'][:])
# add DM assembly
grp_DM = lum_file.create_group("DM")
grp_DM.create_dataset('DM_time', data=SFH_file['DM/DM_time'][:])
grp_DM.create_dataset('DM_z', data=SFH_file['DM/DM_z'][:])
grp_DM.create_dataset('DM_M', data=SFH_file['DM/DM_M'][:])
grp_DM.create_dataset('DM_Mt', data=SFH_file['DM/DM_Mt'][:])
# add luminosities
grp_lum = lum_file.create_group("luminosities")

# close other (SFH) hdf5 file
SFH_file.close()


# iterate over all files, create dictionary with luminosity estimates

dict_lum_data = {}

for ii_file in range(number_of_bins):
    SP_file = h5py.File(path_SP_cat + '/' + filename_SFH_file[:-5] + '/' + filename_SP_file[:-5] + '_' + str(ii_file) + '.hdf5', 'r')
    if (ii_file == 0):
        for ii_key in SP_file['luminosities'].keys():
            if 'luminosity' in ii_key:
                dict_lum_data[ii_key] = SP_file['luminosities'][ii_key][:]
    else:
        for ii_key in SP_file['luminosities'].keys():
            if 'luminosity' in ii_key:
                dict_lum_data[ii_key] = np.vstack([dict_lum_data[ii_key], SP_file['luminosities'][ii_key][:]])
    SP_file.close()


# add this dictionary to the hdf5 file

ii_file = 0
SP_file = h5py.File(path_SP_cat + '/' + filename_SFH_file[:-5] + '/' + filename_SP_file[:-5] + '_' + str(ii_file) + '.hdf5', 'r')

# copy attibutes
for ii_key_attr in SP_file['luminosities'].attrs.keys():
    grp_lum.attrs[ii_key_attr] = SP_file['luminosities'].attrs[ii_key_attr]
# copy content
for ii_key in SP_file['luminosities'].keys():
    if 'luminosity' in ii_key:
        subgrp_lum = grp_lum.create_dataset(ii_key, data=dict_lum_data[ii_key])
        for ii_key_attr in SP_file['luminosities'][ii_key].attrs.keys():
            subgrp_lum.attrs[ii_key_attr] = SP_file['luminosities'][ii_key].attrs[ii_key_attr]

SP_file.close()


lum_file.close()


