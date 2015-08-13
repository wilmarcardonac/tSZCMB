#!/bin/sh
#SBATCH --cpus-per-task=16
#SBATCH --job-name=ylMz
#SBATCH --ntasks=1
#SBATCH --time=0-00:15:00
#SBATCH --mail-user=wilmar.cardona@unige.ch
#SBATCH --mail-type=ALL
#SBATCH --partition=debug
#SBATCH --clusters=baobab
#SBATCH --output=slurm-%J.out

srun ./computeylMz