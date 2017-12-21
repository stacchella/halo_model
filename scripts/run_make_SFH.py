'''
Sandro Tacchella
December 20, 2017
'''

# import modules

import numpy as np
import os
import h5py

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

DM_accretion_history_filename = 'MergerHistory_COLOR_CDM_z3.96.hdf5'
efficency_filename = 'Calibration/epsilon.npy'
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

# loop over all halos

for idx_h in range(len(M_table_in))[::1000]:
    print 'progress (%): ', round(100.0*idx_h/len(M_table_in), 3)
    time_list, SFR_list = make_SFH.construct_SFH(Mt_table_in[idx_h], t_snapshots, look_back=round(t_snapshots[0]), dt=0.1, SFH_type=SFH_type_option, epsilon_fct=epsilon_efficency_fct)
    if (idx_h == 0):
        SFH_table_SFR = SFR_list
    else:
        SFH_table_SFR = np.vstack([SFH_table_SFR, SFR_list])


# save SFH

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


