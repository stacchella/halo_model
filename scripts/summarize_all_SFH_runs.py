'''
Sandro Tacchella
January 2, 2018   : iniate
run this script after run_make_SFH.py in order to
combine all (numpy) files to one (hdf5) file.

'''

# import modules

import numpy as np
import argparse
import os
import h5py

import read_in_halo_cat

from astropy.cosmology import WMAP7 as cosmo


# define paths

path_main = os.environ['WDIR_halo_model']
path_DM_cat = path_main + 'catalogs/DM/'
path_SFH_cat = path_main + 'catalogs/SFH/'


# get number of bins

parser = argparse.ArgumentParser()
parser.add_argument("--number_of_bins", type=int, help="number of cores")
parser.add_argument("--filename_SFH", type=str, help="filename of SFH file")
parser.add_argument("--filename_DM", type=str, help="filename of DM file")
args = parser.parse_args()

number_of_bins = args.number_of_bins


# iterate over all files

counter = 0

for ii in range(number_of_bins):
    file_name = path_SFH_cat + '/' + args.filename_SFH[:-5] + '/' + args.filename_SFH[:-5] + '_' + str(int(float(ii))) + '.npy'
    if (counter == 0):
        SFH_table_SFR = np.load(file_name)
        counter += 1
    else:
        SFH_table_SFR = np.vstack([SFH_table_SFR, np.load(file_name)])
        counter += 1
    print 'progress (%): ', round(100.0*counter/number_of_bins, 3)


time_list = np.load(path_SFH_cat + '/' + args.filename_SFH[:-5] + '/' + args.filename_SFH[:-5] + '_t_' + str(int(float(ii))) + '.npy')


# get dark matter accretion history

z_table_in, M_table_in, Mt_table_in = read_in_halo_cat.read_in_halo_cat(path_DM_cat + args.filename_DM, cosmo)

t_snapshots = 10**3*cosmo.age(z_table_in).value  # in Myr


# save combined array as hdf5 file

try:
    os.remove(path_SFH_cat + args.filename_SFH)
except OSError:
    pass


f = h5py.File(path_SFH_cat + args.filename_SFH, 'w')

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



