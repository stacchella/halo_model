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
path_DM_cat = path_main + 'catalogs/DM/'
path_SFH_cat = path_main + 'catalogs/SFH/'
path_SP_cat = path_main + 'catalogs/SP/'


# set parameters

# SFH_type_option = 'constant'  # 'constant' or 'random'

# z=4
# filename_SFH_file = 'SFH_z4_' + SFH_type_option + '.hdf5'
# filename_SP_file = 'SFH_z4_' + SFH_type_option + '_with_L.hdf5'
# filename_SFH_file = 'SFH_z4_' + SFH_type_option + '_calibration.hdf5'
# filename_SP_file = 'SFH_z4_' + SFH_type_option + '_calibration_with_L.hdf5'
# z=6
# filename_SFH_file = 'SFH_z6_' + SFH_type_option + '.hdf5'
# filename_SP_file = 'SFH_z6_' + SFH_type_option + '_with_L.hdf5'
# z=8
# filename_SFH_file = 'SFH_z8_' + SFH_type_option + '.hdf5'
# filename_SP_file = 'SFH_z8_' + SFH_type_option + '_with_L.hdf5'
# z=10
# filename_SFH_file = 'SFH_z10_' + SFH_type_option + '.hdf5'
# filename_SP_file = 'SFH_z10_' + SFH_type_option + '_with_L.hdf5'


# get number of bins

parser = argparse.ArgumentParser()
parser.add_argument("--number_of_bins", type=int, help="number of cores")
parser.add_argument("--filename_SFH", type=str, help="filename of SFH file")
parser.add_argument("--filename_SP", type=str, help="filename of SP file")
args = parser.parse_args()

number_of_bins = args.number_of_bins


# read in SFH

SFH_file = h5py.File(path_SFH_cat + args.filename_SFH, 'r')
SFH_time = SFH_file['SFH/SFH_time'][:]
SFH_SFR = SFH_file['SFH/SFH_SFR'][:]


# set up new hdf5 file for saving luminosities optained from SFH

try:
    os.remove(path_SP_cat + args.filename_SP)
except OSError:
    pass


lum_file = h5py.File(path_SP_cat + args.filename_SP, 'w')
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
grp_lum = lum_file.create_group("SP")
grp_EmL_lum = grp_lum.create_group("EmL")
grp_EmL_lum.attrs['EL_info'] = 'L_Lya', 'L_HeII', 'L_OIII_L1', 'L_OIII_L2', 'L_CIII_1', 'L_CIII_2', 'L_CIV', 'L_OII', 'L_Hb', 'L_OIII', 'L_Ha', 'L_NII', 'L_SII_1', 'L_SII_2'
grp_EmL_lum.attrs['EL_wavelength'] = np.array([1.215670e+03, 1.640420e+03, 1.661240e+03, 1.666150e+03, 1.906680e+03, 1.908730e+03, 1.908730e+03, 3.727100e+03, 4.862710e+03, 5.008240e+03, 6.564600e+03, 6.585270e+03, 6.718290e+03, 6.732670e+03])
grp_FilL_lum = grp_lum.create_group("FilL")
grp_FilL_lum.attrs['FL_info'] = 'i1500', 'i2300', 'i2800', 'v', 'u', '2mass_j'
grp_spec_lum = grp_lum.create_group("spec")



grp_lum = lum_file.create_group("luminosities")

# close other (SFH) hdf5 file
SFH_file.close()


# iterate over all files, create dictionary with luminosity estimates

dict_Ms_data = {}
dict_FL_data = {}
dict_EL_data = {}
dict_spec_data = {}
stellar_mass_data = []
wavelength_data = []

dict_lum_data = {}

for ii_file in range(number_of_bins):
    SP_file = h5py.File(path_SP_cat + '/' + args.filename_SFH[:-5] + '/' + args.filename_SP[:-5] + '_' + str(ii_file) + '.hdf5', 'r')
    if (ii_file == 0):
        for ii_key in SP_file['SP/FilL'].keys():
            if 'luminosity' in ii_key:
                dict_FL_data[ii_key] = SP_file['SP/FilL'][ii_key][:]
    else:
        for ii_key in SP_file['SP/FilL'].keys():
            if 'luminosity' in ii_key:
                dict_FL_data[ii_key] = np.vstack([dict_FL_data[ii_key], SP_file['SP/FilL'][ii_key][:]])
    if (ii_file == 0):
        for ii_key in SP_file['SP/EmL'].keys():
            if 'luminosity' in ii_key:
                dict_EL_data[ii_key] = SP_file['SP/EmL'][ii_key][:]
    else:
        for ii_key in SP_file['SP/EmL'].keys():
            if 'luminosity' in ii_key:
                dict_EL_data[ii_key] = np.vstack([dict_EL_data[ii_key], SP_file['SP/EmL'][ii_key][:]])
    if (ii_file == 0):
        for ii_key in SP_file['SP/spec'].keys():
            if 'luminosity' in ii_key:
                dict_spec_data[ii_key] = SP_file['SP/spec'][ii_key][:]
    else:
        for ii_key in SP_file['SP/spec'].keys():
            if 'luminosity' in ii_key:
                dict_spec_data[ii_key] = np.vstack([dict_spec_data[ii_key], SP_file['SP/spec'][ii_key][:]])



    SP_file.close()


# add this dictionary to the hdf5 file

ii_file = 0
SP_file = h5py.File(path_SP_cat + '/' + args.filename_SFH[:-5] + '/' + args.filename_SP[:-5] + '_' + str(ii_file) + '.hdf5', 'r')

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


