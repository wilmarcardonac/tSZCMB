#!/bin/sh
#SBATCH --cpus-per-task=16
#SBATCH --job-name=ylMz
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=wilmar.cardona@unige.ch
#SBATCH --mail-type=ALL
#SBATCH --partition=dpt
#SBATCH --clusters=baobab
#SBATCH --output=slurm-%J.out

srun ./computeylMz