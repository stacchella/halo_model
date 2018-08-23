#!/bin/bash
### Name of the job
### Requested number of cores
#SBATCH -n 1
### Requested number of nodes
#SBATCH -N 1
### Requested computing time in minutes
#SBATCH -t 10080
### Partition or queue name
#SBATCH -p conroy
### memory per cpu, in MB
#SBATCH --mem-per-cpu=4000
### Job name
#SBATCH -J 'SFH_z4_calibration'
### output and error logs
#SBATCH -o SFH_z4_calibration_%a.out
#SBATCH -e SFH_z4_calibration_%a.err
### mail
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.tacchella@cfa.harvard.edu
source activate pro
srun -n 1 python /n/eisenstein_lab/Users/stacchella/halo_model/scripts/run_make_SFH.py \
--number_of_bins=200 \
--idx_halo_key="${SLURM_ARRAY_TASK_ID}" \
--SFH_type="constant" \
--filename_SFH="SFH_z4_calibration.hdf5" \
--redshift=4 \
--calibration_run="True" \
--filename_efficiency="None" \
