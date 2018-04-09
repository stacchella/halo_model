# halo_model

This is a simple model that links the dark matter halo growth to the star-formation rate of its galaxy. 

/scripts

1  run run_make_SFH.py to generate SFH from halo catalog (make sure directories are created)
   sbatch --array=1-100 submission_script_make_SFH.sh, with XX given by number_of_bins

2  summarize the run with summarize_all_SFH_runs.py
   run with python pro in cluster environment with large memory (start_interactive_session_largeMEM, source activate pro)
   python summarize_all_SFH_runs.py --number_of_bins 100 --filename_SFH SFH_z10_constant.hdf5 --filename_DM MergerHistory_COLOR_CDM_z10.00.hdf5

3  run run_derive_SP.py to compute the luminosities for the given SFHs
   sbatch --array=1-200 submission_script_derive_SP.sh, with XX given by number_of_bins

4  summarize the run with summarize_all_SP_runs.py
   python summarize_all_SP_runs.py --number_of_bins 200 --filename_SFH SFH_z10_constant.hdf5 --filename_SP SFH_z10_constant_with_L.hdf5


logbook of runs:


- calibration run:  run SFR=1 in order to get a first guess on efficency
    sbatch --array=1-100 submission_script_calibration_SFH.sh
    python summarize_all_SFH_runs.py --number_of_bins 100 --filename_SFH SFH_z4_calibration.hdf5 --redshift 4
    sbatch --array=1-200 submission_script_calibration_SP.sh
    python summarize_all_SP_runs.py --number_of_bins 200 --filename_SFH SFH_z4_calibration.hdf5 --filename_SP snapshot_z4_calibration.hdf5

- tuning fiducial metallicity of 0.02 Zsun

    sbatch --array=1-100 submission_script_tuning_SFH.sh (update file names)
    python summarize_all_SFH_runs.py --number_of_bins 100 --filename_SFH SFH_z4_Z_fid_0X.hdf5 --redshift 4
    sbatch --array=1-200 submission_script_tuning_SP.sh (update file names)
    python summarize_all_SP_runs.py --number_of_bins 200 --filename_SFH SFH_z4_Z_fid_0X.hdf5 --filename_SP snapshot_z4_Z_fid_0X.hdf5

    scp -r stacchella@odyssey.rc.fas.harvard.edu://n/eisenstein_lab/Users/stacchella/halo_model/catalogs/SP/*.hdf5 /Volumes/Tacchella/Work/Postdoc/Halo_Model/snapshots/


- tunning run:		change high end of efficiency in order to match UV LF
					we run this 3 times, each time with different metallicity (Zsun, 0.2Zsun, 0.02Zsun) but with same IMF (Salpeter IMF)
					0.02 Zsun: run_tunning_Z_low_XX.hdf5
					0.2 Zsun:  run_tunning_Z_fid_XX.hdf5
					1.0 Zsun:  run_tunning_Z_hig_XX.hdf5
					where XX can be the subruns for the different high end slopes (actually tunning)
					if necessary and of interest, vary dust prescription in the same way as metallicity.

- actual runs:		use the three different efficiencies from the varying Z runs
					run_Z_Y_zX.hdf5  with Y = {low, fid, hig}, X = {4, 6, 8, 10, ...}

python summarize_all_SFH_runs.py --number_of_bins 100 --filename_SFH SFH_z4_Z_fid.hdf5 --filename_DM MergerHistory_COLOR_CDM_z3.96.hdf5

python summarize_all_SP_runs.py --number_of_bins 200 --filename_SFH SFH_z4_Z_fid.hdf5 --filename_SP snapshot_z4_Z_fid.hdf5


