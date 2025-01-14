#!/bin/bash

#SBATCH --job-name urea_1.3
#SBATCH --nodes=1
#SBATCH --tasks=15
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=10:30:00
#SBATCH --account=iontransport
#SBATCH --error=std.err_%j
#SBATCH --output=std.out_%j
#SBATCH --partition shared

# Load CP2K
module purge
module use /nopt/nrel/apps/modules/centos77/modulefiles/
module load mpich
module load intel
module load cp2k

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "begin job.."
echo $PWD

srun cp2k.psmp -i geo_opt.inp -o out_opt.log
wait
srun cp2k.psmp -i main.inp -o out_main.log
