'''
Sandro Tacchella
January 2, 2018   : iniate
run this script after run_make_SFH.py in order to
combine all (numpy) files to one (hdf5) file.

'''

# import modules

import numpy as np
import glob
import os
import h5py

import read_in_halo_cat

from astropy.cosmology import WMAP7 as cosmo


# define paths

path_main = os.environ['WDIR_halo_model']
path_figures = path_main + 'Figures/Tests/'
path_DM_cat = path_main + 'catalogs/DM/'
path_SFH_cat = path_main + 'catalogs/SFH/'


# set parameters

DM_accretion_history_filename = 'MergerHistory_COLOR_CDM_z3.96.hdf5'
filename_SFH_file = 'SFH_z4_random.hdf5'


# get all files

list_files_all = glob.glob(path_SFH_cat + '/' + filename_SFH_file[:-5] + '/' + filename_SFH_file[:-5] + '*.npy')


# iterate over it

counter = 0
time_list = None

for file_name in list_files_all:
    if ('_t_' in file_name) and (time_list is None):
        time_list = np.load(file_name)
    elif (counter == 0) and ('_t_' not in file_name):
        SFH_table_SFR = np.load(file_name)
        counter += 1
    elif ('_t_' not in file_name):
        SFH_table_SFR = np.vstack([SFH_table_SFR, np.load(file_name)])
        counter += 1
    print 'progress (%): ', round(100.0*counter/(0.5*len(list_files_all)), 3)


# get dark matter accretion history

z_table_in, M_table_in, Mt_table_in = read_in_halo_cat.read_in_halo_cat(path_DM_cat + DM_accretion_history_filename, cosmo)

t_snapshots = 10**3*cosmo.age(z_table_in).value  # in Myr


# save combined array as hdf5 file

try:
    os.remove(path_SFH_cat + filename_SFH_file)
except OSError:
    pass


f = h5py.File(path_SFH_cat + filename_SFH_file, 'w')

# add SFH
grp_SFH = f.create_group("SFH")
grp_SFH.create_dataset('SFH_time', data=time_list)
grp_SFH.create_dataset('SFH_SFR', data=SFH_table_SFR)

# add DM assembly
grp_DM = f.create_group("DM")
grp_DM.create_dataset('DM_time', data=t_snapshots)
grp_DM.create_dataset('DM_z', data=z_table_in)
grp_DM.create_dataset('DM_M', data=M_table_in)
grp_DM.create_dataset('DM_Mt', data=Mt_table_in)
f.close()


