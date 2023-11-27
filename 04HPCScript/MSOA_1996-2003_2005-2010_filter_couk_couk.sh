#!/bin/bash

#SBATCH --job-name=test_MSOA_1996-2003_2005-2010_filter_couk-couk
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=08:30:00
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

Rscript 1996-2003_2005-2010_couk_couk_filter.R

