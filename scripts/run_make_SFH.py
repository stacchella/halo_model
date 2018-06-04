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
path_SFH_cat = path_main + 'catalogs/SFH/'


# read in command line arguments

parser = argparse.ArgumentParser()
parser.add_argument("--number_of_bins", type=int, help="number of cores")
parser.add_argument("--idx_halo_key", type=int, help="iteration variable")
parser.add_argument("--SFH_type", type=str, help="SFH type: constant or random")
parser.add_argument("--filename_SFH", type=str, help="filename of SFH file")
parser.add_argument("--redshift", type=int, help="redshift")
parser.add_argument("--calibration_run", type=str, help="True or False")
parser.add_argument("--filename_efficiency", type=str, help="filename of efficiency file")
args = parser.parse_args()


run_params = {'number_of_bins': args.number_of_bins,  # this gives number of cores we run on
              'idx_halo_key': args.idx_halo_key,  # iteration variable
              }

# run_params = {'number_of_bins': 1,  # this gives number of cores we run on
#              'idx_halo_key': 1,  # iteration variable
#              }


# get dark matter accretion history, remove contaminated halos

z_table_in, M_table_in, Mt_table_in, is_contam = read_in_halo_cat.read_in_halo_cat(args.redshift, cosmo)


print len(z_table_in)
print len(M_table_in)
print len(Mt_table_in)


# set up efficency function (based on calibration)

if (args.calibration_run == 'True'):
    def epsilon_efficency_fct(Mh_in, size_in=1.0):
        '''
        This function returns an efficency from
        the calibrated distribution for a given halo mass.
        '''
        return(np.zeros(size_in))
else:
    epsilon_efficency_fct = read_in_efficency.read_in_efficency(path_SFH_cat + args.filename_efficiency)


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
    time_list, SFR_list, Mz_list, Z_list = make_SFH.construct_SFH(Mt_table_in[idx_h], t_snapshots, SFH_type=args.SFH_type, epsilon_fct=epsilon_efficency_fct, dt_high_res=0.1, dt_low_res=20.0, time_delay=0.1, time_smoothing=0.02, specific_growth_threshold=1.0, Z0=0.0143*10**-3, R=0.1, y=0.023, lam10=0.4)
    if (counter == 0):
        SFH_table_SFR = SFR_list
        SFH_table_Mz = Mz_list
        SFH_table_Z = Z_list
    else:
        SFH_table_SFR = np.vstack([SFH_table_SFR, SFR_list])
        SFH_table_Mz = np.vstack([SFH_table_Mz, Mz_list])
        SFH_table_Z = np.vstack([SFH_table_Z, Z_list])
    counter += 1


# save SFH as numpy file (in a new directory), later combine all these files

np.save(path_SFH_cat + '/' + args.filename_SFH[:-5] + '/' + args.filename_SFH[:-5] + '_' + str(int(float(args.idx_halo_key))-1) + '.npy', SFH_table_SFR)
np.save(path_SFH_cat + '/' + args.filename_SFH[:-5] + '/' + args.filename_SFH[:-5] + '_Mz_' + str(int(float(args.idx_halo_key))-1) + '.npy', SFH_table_Mz)
np.save(path_SFH_cat + '/' + args.filename_SFH[:-5] + '/' + args.filename_SFH[:-5] + '_Z_' + str(int(float(args.idx_halo_key))-1) + '.npy', SFH_table_Z)
np.save(path_SFH_cat + '/' + args.filename_SFH[:-5] + '/' + args.filename_SFH[:-5] + '_t_' + str(int(float(args.idx_halo_key))-1) + '.npy', time_list)

#np.save(path_SFH_cat + '/' + args.filename_SFH[:-5] + '/' + args.filename_SFH[:-5] + '_0.npy', SFH_table_SFR)
#np.save(path_SFH_cat + '/' + args.filename_SFH[:-5] + '/' + args.filename_SFH[:-5] + '_t_0.npy', time_list)


