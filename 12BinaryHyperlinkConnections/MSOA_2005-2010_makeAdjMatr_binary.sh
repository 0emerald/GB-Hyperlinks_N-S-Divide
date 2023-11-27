#!/bin/bash

#SBATCH --job-name=MSOA_2005-2010_makeAdjMatr_binary
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=07:30:00
#SBATCH --mem=128G

cd "${SLURM_SUBMIT_DIR}"

echo Running on host "$(hostname)"
echo Time is "$(date)"
echo Directory is "$(pwd)"
echo Slurm job ID is "${SLURM_JOBID}"
echo This jobs runs on the following machines:
echo "${SLURM_JOB_NODELIST}"

module add lang/r/4.1.2-gcc
module load lang/r/4.1.2-gcc

export OMP_NUM_THREADS=1

Rscript MSOA_2005-2010_makeAdjMatr_binary.R

