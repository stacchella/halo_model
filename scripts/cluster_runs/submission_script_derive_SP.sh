#!/bin/bash
### Name of the job
### Requested number of cores
#SBATCH -n 1
### Requested number of nodes
#SBATCH -N 1
### Requested computing time in minutes
#SBATCH -t 10080
### Partition or queue name
#SBATCH -p general
### memory per cpu, in MB
#SBATCH --mem-per-cpu=60000
### Job name
#SBATCH -J 'random_SFH'
### output and error logs
#SBATCH -o random_SFH_L_%a.out
#SBATCH -e random_SFH_L_%a.err
### mail
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.tacchella@cfa.harvard.edu
source activate pro
srun -n 1 python /n/regal/eisenstein_lab/stacchella/halo_model/scripts/run_derive_SP.py \
--number_of_bins=200 \
--idx_halo_key="${SLURM_ARRAY_TASK_ID}"

