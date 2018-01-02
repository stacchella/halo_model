'''
Sandro Tacchella
December 20, 2017 : iniate
January 2, 2018   : update parallel
=> sbatch --array=1-21 submission_script_make_SFH.sh

'''

# import modules

import numpy as np
import os

import read_in_halo_cat
import read_in_efficency
import make_SFH

from astropy.cosmology import WMAP7 as cosmo


# define parameters

run_params = {'number_of_bins': 20,  # this gives number of cores we run on
              'idx_halo_key': 0.0,  # iteration variable
              }


# define paths

path_main = os.environ['WDIR_halo_model']
path_figures = path_main + 'Figures/Tests/'
path_DM_cat = path_main + 'catalogs/DM/'
path_SFH_cat = path_main + 'catalogs/SFH/'


# set parameters

DM_accretion_history_filename = 'MergerHistory_COLOR_CDM_z3.96.hdf5'
efficency_filename = 'calibration/epsilon.npy'
SFH_type_option = 'random'
filename_SFH_file = 'SFH_z4_random.hdf5'


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

def get_halo_ids(number_of_bins, idx_halo_key=1.0):
    idx_all_halos = range(len(M_table_in))
    idx_bins_all_halos = np.array_split(idx_all_halos, number_of_bins)
    return(idx_bins_all_halos[int(idx_halo_key-1)])  # -1 since slurm counts from 1 (and not from 0)


idx_halo_considered = get_halo_ids(**run_params)


# loop over all halos

for idx_h in idx_halo_considered[::100]:
    print 'progress (%): ', round(100.0*idx_h/len(M_table_in), 3)
    time_list, SFR_list = make_SFH.construct_SFH(Mt_table_in[idx_h], t_snapshots, look_back=round(t_snapshots[0]), dt=0.1, SFH_type=SFH_type_option, epsilon_fct=epsilon_efficency_fct)
    if (idx_h == 0):
        SFH_table_SFR = SFR_list
    else:
        SFH_table_SFR = np.vstack([SFH_table_SFR, SFR_list])


# save SFH as numpy file, later combine all these files

np.save(path_SFH_cat + filename_SFH_file[:-5] + '_' + str(run_params['idx_halo_key']) + '.npy', SFH_table_SFR)
np.save(path_SFH_cat + filename_SFH_file[:-5] + '_t_' + str(run_params['idx_halo_key']) + '.npy', time_list)



