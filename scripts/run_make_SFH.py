'''
Sandro Tacchella
December 20, 2017 : iniate
January 2, 2018   : update parallel
=> sbatch --array=1-XX submission_script_make_SFH.sh, with XX given by number_of_bins

'''

# import modules

import numpy as np
import os
import argparse

import read_in_halo_cat
import read_in_efficency
import make_SFH

from astropy.cosmology import WMAP7 as cosmo


# define paths

path_main = os.environ['WDIR_halo_model']
path_figures = path_main + 'Figures/Tests/'
path_DM_cat = path_main + 'catalogs/DM/'
path_SFH_cat = path_main + 'catalogs/SFH/'


# set parameters

efficency_filename = 'calibration/epsilon.npy'
SFH_type_option = 'random'
# z=4
DM_accretion_history_filename = 'MergerHistory_COLOR_CDM_z3.96.hdf5'
filename_SFH_file = 'SFH_z4_random.hdf5'
# z=6
#DM_accretion_history_filename = 'MergerHistory_COLOR_CDM_z5.98.hdf5'
#filename_SFH_file = 'SFH_z6_random.hdf5'
# z=8
#DM_accretion_history_filename = 'MergerHistory_COLOR_CDM_z8.10.hdf5'
#filename_SFH_file = 'SFH_z8_random.hdf5'
# z=10
#DM_accretion_history_filename = 'MergerHistory_COLOR_CDM_z10.00.hdf5'
#filename_SFH_file = 'SFH_z10_random.hdf5'


# read in command line arguments

parser = argparse.ArgumentParser()
parser.add_argument("--number_of_bins", type=int, help="number of cores")
parser.add_argument("--idx_halo_key", type=int, help="iteration variable")
args = parser.parse_args()


run_params = {'number_of_bins': args.number_of_bins,  # this gives number of cores we run on
              'idx_halo_key': args.idx_halo_key,  # iteration variable
              }

# run_params = {'number_of_bins': 1,  # this gives number of cores we run on
#              'idx_halo_key': 1,  # iteration variable
#              }


# get dark matter accretion history

z_table_in, M_table_in, Mt_table_in = read_in_halo_cat.read_in_halo_cat(path_DM_cat + DM_accretion_history_filename, cosmo)

print len(z_table_in)
print len(M_table_in)
print len(Mt_table_in)


# set up efficency function (based on calibration)

epsilon_efficency_fct = read_in_efficency.read_in_efficency(path_SFH_cat + efficency_filename)


# get SFH: random burst in last step

t_snapshots = 10**3*cosmo.age(z_table_in).value  # in Myr


# split halo in bins

def get_halo_ids(number_of_bins, idx_halo_key=1.0, **kwargs):
    idx_all_halos = range(len(M_table_in))
    idx_bins_all_halos = np.array_split(idx_all_halos, number_of_bins)
    print idx_bins_all_halos[int(float(idx_halo_key))-1]
    return(idx_bins_all_halos[int(float(idx_halo_key))-1])  # -1 since slurm counts from 1 (and not from 0)


idx_halo_considered = get_halo_ids(**run_params)


# loop over all halos

counter = 0

for idx_h in idx_halo_considered:
    print 'progress (%): ', round(100.0*counter/len(idx_halo_considered), 3)
    if (round(t_snapshots[0]) > 200.0):
        time_list_highres, SFR_list_highres = make_SFH.construct_SFH(Mt_table_in[idx_h], t_snapshots, look_back=200.0, dt=0.1, SFH_type=SFH_type_option, epsilon_fct=epsilon_efficency_fct)
        time_list_lowres, SFR_list_lowres = make_SFH.construct_SFH(Mt_table_in[idx_h], t_snapshots, look_back=round(t_snapshots[0]), dt=20.0, SFH_type='constant', epsilon_fct=epsilon_efficency_fct)
        idx = (np.abs(time_list_lowres-time_list_highres[2])).argmin()
        time_list = np.append(time_list_lowres[:idx], time_list_highres[2:])
        SFR_list = np.append(SFR_list_lowres[:idx], SFR_list_highres[2:])
    else:
        time_list, SFR_list = make_SFH.construct_SFH(Mt_table_in[idx_h], t_snapshots, look_back=round(t_snapshots[0]), dt=0.1, SFH_type=SFH_type_option, epsilon_fct=epsilon_efficency_fct)
    if (counter == 0):
        SFH_table_SFR = SFR_list
    else:
        SFH_table_SFR = np.vstack([SFH_table_SFR, SFR_list])
    counter += 1


# save SFH as numpy file (in a new directory), later combine all these files

np.save(path_SFH_cat + '/' + filename_SFH_file[:-5] + '/' + filename_SFH_file[:-5] + '_' + str(int(float(args.idx_halo_key))-1) + '.npy', SFH_table_SFR)
np.save(path_SFH_cat + '/' + filename_SFH_file[:-5] + '/' + filename_SFH_file[:-5] + '_t_' + str(int(float(args.idx_halo_key))-1) + '.npy', time_list)

#np.save(path_SFH_cat + '/' + filename_SFH_file[:-5] + '/' + filename_SFH_file[:-5] + '_0.npy', SFH_table_SFR)
#np.save(path_SFH_cat + '/' + filename_SFH_file[:-5] + '/' + filename_SFH_file[:-5] + '_t_0.npy', time_list)


