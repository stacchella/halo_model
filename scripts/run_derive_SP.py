'''
Sandro Tacchella
December 21, 2017 : intiate
January 2, 2018   : update parallel


NOTES:

We start with a fiducial model:
- IMF         :   Salpeter (imf_type=0)
- metallicity :   solar metallicity (logzsol=0.0)
- dust        :   average beta relation

Possible variations:
- IMF         :
                * 0: Salpeter (1955)
                * 1: Chabrier (2003)
                * 2: Kroupa (2001)
                * 3: van Dokkum (2008)
                * 4: Dave (2008)
                * 5: tabulated piece-wise power law IMF
- metallicity : simple grid
- dust        : apply scatter; add dust to nebular in some peculiar way ("fudge factor" of 0.44 to 1.00)

Things to keep in mind:
- stocastics effects at very low SFRs

'''

# import modules

import numpy as np
import os
import h5py
import fsps
import itertools

import derive_SP_prop


# define paths

path_main = os.environ['WDIR_halo_model']
path_figures = path_main + 'Figures/Tests/'
path_SFH_cat = path_main + 'catalogs/SFH/'
path_SP_cat = path_main + 'catalogs/SP/'


# set parameters

filename_SFH_file = 'SFH_z4_random.hdf5'
filename_SP_file = 'SFH_z4_random_with_L.hdf5'
idx_every_other = 2  # fundge factor so that SFH setting works


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
grp_SFH.create_dataset('SFH_SFR', data=SFH_file['SFH/SFH_time'][:])
# add DM assembly
grp_DM = lum_file.create_group("DM")
grp_DM.create_dataset('DM_time', data=SFH_file['DM/DM_time'][:])
grp_DM.create_dataset('DM_z', data=SFH_file['DM/DM_z'][:])
grp_DM.create_dataset('DM_M', data=SFH_file['DM/DM_M'][:])
grp_DM.create_dataset('DM_Mt', data=SFH_file['DM/DM_Mt'][:])
# add luminosities
grp_lum = lum_file.create_group("luminosities")
grp_lum.attrs['luminosity_info'] = 'L_Ha', 'L_Hb', 'L_1500', 'L_2300', 'L_2800'

# close other (SFH) hdf5 file
SFH_file.close()


# set up the grid

logzsol_grid = [0.0, -0.5, -1.0, -1.5, -2.0]
IMF_grid = [0, 1, 3]

iterables = [logzsol_grid, IMF_grid]

counter = 0.0

for t in itertools.product(*iterables):
    if (counter == 0.0):
        all_combinations = t
    else:
        all_combinations = np.vstack([all_combinations, t])
    counter += 1.0


dict_names = ['logzsol', 'IMF']
dict_all_combinations = []

for ii in range(all_combinations.shape[0]):
    ii_dict = {}
    for jj in range(len(dict_names)):
        ii_dict[dict_names[jj]] = all_combinations[ii][jj]
        if (dict_names[jj] == 'IMF'):
            if (all_combinations[ii][jj] == 0.0):
                ii_dict['IMF_name'] = 'Salpeter'
            elif (all_combinations[ii][jj] == 1.0):
                ii_dict['IMF_name'] = 'Chabrier'
            elif (all_combinations[ii][jj] == 3.0):
                ii_dict['IMF_name'] = 'van Dokkum'
    if (ii == 0.0):
        dict_all_combinations = ii_dict
    else:
        dict_all_combinations = np.vstack([dict_all_combinations, ii_dict])


print 'this is the grid for the SP:'
print dict_all_combinations


# iterate over all models and compute luminosities

for ii_model in range(len(dict_all_combinations)):
    model_dict = dict_all_combinations[ii_model][0]
    print '================================================='
    print 'working on model with: ', model_dict['logzsol']
    print '                   log Z/Zsun = ', model_dict['logzsol']
    print '                   IMF        = ', model_dict['IMF_name']
    # set standard SP model
    # here we set metallicity and IMF
    sp_now = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1, imf_type=model_dict['IMF'], add_neb_emission=True, sfh=3, logzsol=model_dict['logzsol'], dust_type=2, dust2=0.0)
    # iterate over all SFHs
    for ii in range(SFH_SFR.shape[0]):
        print 'progress (%): ', round(100.0*ii/SFH_SFR.shape[0], 3)
        lum_vec = derive_SP_prop.get_luminosities_for_SFH(sp_now, [10**-3*SFH_time[::idx_every_other], SFH_SFR[ii][::idx_every_other]], [10**-3*SFH_time[::idx_every_other][-1]])
        lum_vec = np.array(lum_vec).T
        if (ii == 0):
            lum_mat = lum_vec
        else:
            lum_mat = np.vstack([lum_mat, lum_vec])
    # add SFH
    subgrp_lum = grp_lum.create_dataset('luminosity_' + str(ii_model), data=lum_mat)
    for ii_key in model_dict.keys():
        subgrp_lum.attrs[ii_key] = model_dict[ii_key]


lum_file.close()





