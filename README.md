# halo_model

This is a simple model that links the dark matter halo growth to the star-formation rate of its galaxy. 

/scripts

1  run run_make_SFH.py to generate SFH from halo catalog
   sbatch --array=1-XX submission_script_make_SFH.sh, with XX given by number_of_bins

2  summarize the run with summarize_all_SFH_runs.py
   run with python pro in cluster environment with large memory (start_interactive_session_largeMEM, source activate pro)
   python summarize_all_SFH_runs.py --number_of_bins XX

3  run run_derive_SP.py to compute the luminosities for the given SFHs
   sbatch --array=1-XX submission_script_derive_SP.sh, with XX given by number_of_bins

4  summarize the run with summarize_all_SP_runs.py
   python summarize_all_SP_runs.py --number_of_bins XX

