#!/bin/bash
### Name of the job
### Requested number of cores
#SBATCH -n 1
### Requested number of nodes
#SBATCH -N 1
### Requested computing time in minutes
#SBATCH -t 10080
### Partition or queue name
#SBATCH -p itc_cluster
### memory per cpu, in MB
#SBATCH --mem-per-cpu=4000
### Job name
#SBATCH -J 'SFH_z4_Z_fid_06'
### output and error logs
#SBATCH -o SFH_z4_Z_fid_06_%a.out
#SBATCH -e SFH_z4_Z_fid_06_%a.err
### mail
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.tacchella@cfa.harvard.edu
source activate pro
srun -n 1 python /n/eisenstein_lab/Users/stacchella/halo_model/scripts/run_make_SFH.py \
--number_of_bins=200 \
--idx_halo_key="${SLURM_ARRAY_TASK_ID}" \
--SFH_type="constant" \
--filename_SFH="SFH_z4_Z_fid_06.hdf5" \
--redshift=4 \
--calibration_run="False" \
--filename_efficiency="calibration/epsilon_constant_param_0.6.npy" \
