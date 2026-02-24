#!/bin/bash

#SBATCH --job-name Urea_1.0
#SBATCH --nodes=4
#SBATCH --tasks-per-node=16
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4G
#SBATCH --time=36:30:00
#SBATCH --account=iontransport
#SBATCH --error=std.err_%j
#SBATCH --output=std.out_%j
#SBATCH --partition medmem

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
srun cp2k.psmp -i nvt_main.inp -o out_nvt.log
wait
srun cp2k.psmp -i nvt_highT.inp -o out_nvt_highT.log
wait
srun cp2k.psmp -i npt_main.inp -o out_npt.log
wait
srun cp2k.psmp -i nvt_prod.inp -o out_nvt_prod.log
