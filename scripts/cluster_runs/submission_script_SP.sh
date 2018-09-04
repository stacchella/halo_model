#!/bin/bash
### Name of the job
### Requested number of cores
#SBATCH -n 1
### Requested number of nodes
#SBATCH -N 1
### Requested computing time in minutes
#SBATCH -t 10080
### Partition or queue name
#SBATCH -p conroy,itc_cluster,hernquist
### memory per cpu, in MB
#SBATCH --mem-per-cpu=4000
### Job name
#SBATCH -J 'SP_z8_Z_fid_02'
### output and error logs
#SBATCH -o SP_z8_Z_fid_02_%a.out
#SBATCH -e SP_z8_Z_fid_02_%a.err
### mail
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.tacchella@cfa.harvard.edu
source activate pro
srun -n 1 python /n/eisenstein_lab/Users/stacchella/halo_model/scripts/run_derive_SP.py \
--number_of_bins=200 \
--idx_halo_key="${SLURM_ARRAY_TASK_ID}" \
--filename_SFH="SFH_z8_Z_fid_02.hdf5" \
--filename_SP="snapshot_z8_Z_fid_02.hdf5" \
--calibration_run="False" \
